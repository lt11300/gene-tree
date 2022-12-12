# gene-tree
import os
import contextlib
import itertools
import collections
import re
import copy
import random
import warnings
from io import StringIO

#基础的树元素

def _level_traverse(root, get_children):
    Q = collections.deque([root])
    while Q:
        v = Q.popleft()
        yield v
        Q.extend(get_children(v))


def _preorder_traverse(root, get_children):

    def dfs(elem):
        yield elem
        for v in get_children(elem):
            yield from dfs(v)

    yield from dfs(root)


def _postorder_traverse(root, get_children):
   
    def dfs(elem):
        for v in get_children(elem):
            yield from dfs(v)
        yield elem

    yield from dfs(root)
    
def _sorted_attrs(elem):
   
    singles = []
    lists = []
    # Sort attributes for consistent results
    for attrname, child in sorted(elem.__dict__.items(), key=lambda kv: kv[0]):
        if child is None:
            continue
        if isinstance(child, list):
            lists.extend(child)
        else:
            singles.append(child)
    return (x for x in singles + lists if isinstance(x, TreeElement))

class TreeElement:
    def __repr__(self):
        def pair_as_kwarg_string(key, val):
            if isinstance(val, str):
                val = val[:57] + "..." if len(val) > 60 else val
                return "%s='%s'" % (key, val)
            return "%s=%s" % (key, val)

        return "%s(%s)" % (
            self.__class__.__name__,
            ", ".join(
                pair_as_kwarg_string(key, val)
                for key, val in sorted(self.__dict__.items())
                if val is not None and type(val) in (str, int, float, bool, str)
            ),
        )

    __str__ = __repr__
    
class TreeMixin:
    # Traversal methods
    def _filter_search(self, filter_func, order, follow_attrs):
        order_opts = {
            "preorder": _preorder_traverse,
            "postorder": _postorder_traverse,
            "level": _level_traverse,
        }
        try:
            order_func = order_opts[order]
        except KeyError:
            raise ValueError(
                "Invalid order '%s'; must be one of: %s" % (order, tuple(order_opts))
            ) from None

        if follow_attrs:
            get_children = _sorted_attrs
            root = self
        else:
            get_children = lambda elem: elem.clades  # noqa: E731
            root = self.root
        return filter(filter_func, order_func(root, get_children))

    def find_any(self, *args, **kwargs):
        hits = self.find_elements(*args, **kwargs)
        try:
            return next(hits)
        except StopIteration:
            return None

    def find_elements(self, target=None, terminal=None, order="preorder", **kwargs):
        if terminal is not None:
            kwargs["terminal"] = terminal
        is_matching_elem = _combine_matchers(target, kwargs, False)
        return self._filter_search(is_matching_elem, order, True)

    def find_clades(self, target=None, terminal=None, order="preorder", **kwargs):

        def match_attrs(elem):
            orig_clades = elem.__dict__.pop("clades")
            found = elem.find_any(target, **kwargs)
            elem.clades = orig_clades
            return found is not None

        if  terminal is None:
            is_matching_elem = match_attrs
        else:
            def is_matching_elem(elem):
                return (elem.is_terminal() == terminal) and match_attrs(elem)

        return self._filter_search(is_matching_elem, order, False)
    
    def get_nonterminals(self, order="preorder"):
        return list(self.find_clades(terminal=False, order=order))

    def get_terminals(self, order="preorder"):
        return list(self.find_clades(terminal=True, order=order))
#树和分支
class Tree(TreeElement, TreeMixin):
    def __init__(self, root=None, rooted=True, id=None, name=None):
        self.root = root or Clade()
        self.rooted = rooted
        self.id = id
        self.name = name

    @classmethod
    def from_clade(cls, clade, **kwargs):
        root = copy.deepcopy(clade)
        return cls(root, **kwargs)

    @classmethod
    def randomized(cls, taxa, branch_length=1.0, branch_stdev=None):
        if isinstance(taxa, int):
            taxa = ["taxon%s" % (i + 1) for i in range(taxa)]
        elif hasattr(taxa, "__iter__"):
            taxa = list(taxa)
        else:
            raise TypeError(
                "taxa argument must be integer (# taxa) or iterable of taxon names."
            )
        rtree = cls()
        terminals = [rtree.root]
        while len(terminals) < len(taxa):
            newsplit = random.choice(terminals)
            newsplit.split(branch_length=branch_length)
            newterms = newsplit.clades
            if branch_stdev:
                # Add some noise to the branch lengths
                for nt in newterms:
                    nt.branch_length = max(0, random.gauss(branch_length, branch_stdev))
            terminals.remove(newsplit)
            terminals.extend(newterms)
        # Distribute taxon labels randomly
        random.shuffle(taxa)
        for node, name in zip(terminals, taxa):
            node.name = name
        return rtree

    @property
    def clade(self):
        return self.root

    def is_terminal(self):
        """Check if the root of this tree is terminal."""
        return not self.root.clades

    # Convention from SeqRecord and Alignment classes

    def __format__(self, format_spec):
        """Serialize the tree as a string in the specified file format.

        This method supports Python's ``format`` built-in function.

        :param format_spec: a lower-case string supported by ``Bio.Phylo.write``
            as an output file format.

        """
        if format_spec:
            from io import StringIO
            from Bio.Phylo import _io

            handle = StringIO()
            _io.write([self], handle, format_spec)
            return handle.getvalue()
        else:
            # Follow python convention and default to using __str__
            return str(self)

    def format(self, fmt=None, format=None):
        """Serialize the tree as a string in the specified file format.

        :param fmt: a lower-case string supported by ``Bio.Phylo.write``
            as an output file format.

        """
        if format is not None:
            if fmt is not None:
                raise ValueError("The ``format`` argument has been renamed to ``fmt``.")
            warnings.warn(
                "The ``format`` argument has been renamed to ``fmt``.",
                BiopythonDeprecationWarning,
            )
            fmt = format
        return self.__format__(fmt)

    # Pretty-printer for the entire tree hierarchy

    def __str__(self):
        """Return a string representation of the entire tree.

        Serialize each sub-clade recursively using ``repr`` to create a summary
        of the object structure.
        """
        TAB = "    "
        textlines = []

        def print_tree(obj, indent):
            """Recursively serialize sub-elements.

            This closes over textlines and modifies it in-place.
            """
            if isinstance(obj, (Tree, Clade)):
                # Avoid infinite recursion or special formatting from str()
                objstr = repr(obj)
            else:
                objstr = str(obj)
            textlines.append(TAB * indent + objstr)
            indent += 1
            for attr in obj.__dict__:
                child = getattr(obj, attr)
                if isinstance(child, TreeElement):
                    print_tree(child, indent)
                elif isinstance(child, list):
                    for elem in child:
                        if isinstance(elem, TreeElement):
                            print_tree(elem, indent)

        print_tree(self, 0)
        return "\n".join(textlines)
    
class Clade(TreeElement, TreeMixin):

    def __init__(
        self,
        branch_length=None,
        name=None,
        clades=None,
        confidence=None,
        color=None,
        width=None,
    ):
        """Define parameters for the Clade tree."""
        self.branch_length = branch_length
        self.name = name
        self.clades = clades or []
        self.confidence = confidence
        self.color = color
        self.width = width

    @property
    def root(self):
        """Allow TreeMixin methods to traverse clades properly."""
        return self

    def is_terminal(self):
        """Check if this is a terminal (leaf) node."""
        return not self.clades

    # Sequence-type behavior methods

    def __getitem__(self, index):
        """Get clades by index (integer or slice)."""
        if isinstance(index, (int, slice)):
            return self.clades[index]
        ref = self
        for idx in index:
            ref = ref[idx]
        return ref

    def __iter__(self):
        """Iterate through this tree's direct descendent clades (sub-trees)."""
        return iter(self.clades)

    def __len__(self):
        """Return the number of clades directy under the root."""
        return len(self.clades)

    def __bool__(self):
        """Boolean value of an instance of this class (True).

        NB: If this method is not defined, but ``__len__``  is, then the object
        is considered true if the result of ``__len__()`` is nonzero. We want
        Clade instances to always be considered True.
        """
        return True

    def __str__(self):
        """Return name of the class instance."""
        if self.name:
            return self.name[:37] + "..." if len(self.name) > 40 else self.name
        return self.__class__.__name__

#newick格式
class Tree1(Tree):

    def __init__(self, root=None, rooted=False, id=None, name=None, weight=1.0):
        Tree.__init__(
            self, root=root or Clade(), rooted=rooted, id=id, name=name
        )
        self.weight = weight

class Clade1(Clade):
    def __init__(
        self, branch_length=None, name=None, clades=None, confidence=None, comment=None
    ):

        Clade.__init__(
            self,
            branch_length=branch_length,
            name=name,
            clades=clades,
            confidence=confidence,
        )
        self.comment = comment
#读树
tokens = [
    (r"\(", "open parens"),
    (r"\)", "close parens"),
    (r"[^\s\(\)\[\]\'\:\;\,]+", "unquoted node label"),
    (r"\:\ ?[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?", "edge length"),
    (r"\,", "comma"),
    (r"\[(\\.|[^\]])*\]", "comment"),
    (r"\'(\\.|[^\'])*\'", "quoted node label"),
    (r"\;", "semicolon"),
    (r"\n", "newline"),
]
tokenizer = re.compile("(%s)" % "|".join(token[0] for token in tokens))
token_dict = {name: re.compile(token) for token, name in tokens}

def parse(handle, **kwargs):
    return Parser(handle).parse(**kwargs)

def _parse_confidence(text):
    if text.isdigit():
        return int(text)
        # NB: Could make this more consistent by treating as a percentage
        # return int(text) / 100.
    try:
        return float(text)
        # NB: This should be in [0.0, 1.0], but who knows what people will do
        # assert 0 <= current_clade.confidence <= 1
    except ValueError:
        return None


def _format_comment(text):
    return "[%s]" % (text.replace("[", "\\[").replace("]", "\\]"))


def _get_comment(clade):
    try:
        comment = clade.comment
    except AttributeError:
        pass
    else:
        if comment:
            return _format_comment(str(comment))
    return ""

class Parser:
    """Parse a Newick tree given a file handle.

    Based on the parser in ``Bio.Nexus.Trees``.
    """

    def __init__(self, handle):
        """Initialize file handle for the Newick Tree."""
        if handle.read(0) != "":
            raise ValueError("Newick files must be opened in text mode") from None
        self.handle = handle

    @classmethod
    def from_string(cls, treetext):
        """Instantiate the Newick Tree class from the given string."""
        handle = StringIO(treetext)
        return cls(handle)

    def parse(
        self, values_are_confidence=False, comments_are_confidence=False, rooted=False
    ):
        """Parse the text stream this object was initialized with."""
        self.values_are_confidence = values_are_confidence
        self.comments_are_confidence = comments_are_confidence
        self.rooted = rooted
        buf = ""
        for line in self.handle:
            buf += line.rstrip()
            if buf.endswith(";"):
                yield self._parse_tree(buf)
                buf = ""
        if buf:
            # Last tree is missing a terminal ';' character -- that's OK
            yield self._parse_tree(buf)

    def _parse_tree(self, text):
        """Parse the text representation into an Tree object (PRIVATE)."""
        tokens = re.finditer(tokenizer, text.strip())

        new_clade = self.new_clade
        root_clade = new_clade()

        current_clade = root_clade
        entering_branch_length = False

        lp_count = 0
        rp_count = 0
        for match in tokens:
            token = match.group()

            if token.startswith("'"):
                # quoted label; add characters to clade name
                current_clade.name = token[1:-1]

            elif token.startswith("["):
                # comment
                current_clade.comment = token[1:-1]
                if self.comments_are_confidence:
                    # Try to use this comment as a numeric support value
                    current_clade.confidence = _parse_confidence(current_clade.comment)

            elif token == "(":
                # start a new clade, which is a child of the current clade
                current_clade = new_clade(current_clade)
                entering_branch_length = False
                lp_count += 1

            elif token == ",":
                # if the current clade is the root, then the external parentheses
                # are missing and a new root should be created
                if current_clade is root_clade:
                    root_clade = new_clade()
                    current_clade.parent = root_clade
                # start a new child clade at the same level as the current clade
                parent = self.process_clade(current_clade)
                current_clade = new_clade(parent)
                entering_branch_length = False

            elif token == ")":
                # done adding children for this parent clade
                parent = self.process_clade(current_clade)
                if not parent:
                    raise NewickError("Parenthesis mismatch.")
                current_clade = parent
                entering_branch_length = False
                rp_count += 1

            elif token == ";":
                break

            elif token.startswith(":"):
                # branch length or confidence
                value = float(token[1:])
                if self.values_are_confidence:
                    current_clade.confidence = value
                else:
                    current_clade.branch_length = value

            elif token == "\n":
                pass

            else:
                # unquoted node label
                current_clade.name = token

        if not lp_count == rp_count:
            raise NewickError("Number of open/close parentheses do not match.")

        # if ; token broke out of for loop, there should be no remaining tokens
        try:
            next_token = next(tokens)
            raise NewickError(
                "Text after semicolon in Newick tree: %s" % next_token.group()
            )
        except StopIteration:
            pass

        self.process_clade(current_clade)
        self.process_clade(root_clade)
        return Tree1(root=root_clade, rooted=self.rooted)

    def new_clade(self, parent=None):
        """Return new Newick.Clade, optionally with temporary reference to parent."""
        clade = Clade1()
        if parent:
            clade.parent = parent
        return clade

    def process_clade(self, clade):
        """Remove node's parent and return it. Final processing of parsed clade."""
        if (
            (clade.name)
            and not (self.values_are_confidence or self.comments_are_confidence)
            and (clade.confidence is None)
            and (clade.clades)
        ):
            clade.confidence = _parse_confidence(clade.name)
            if clade.confidence is not None:
                clade.name = None

        try:
            parent = clade.parent
        except AttributeError:
            pass
        else:
            parent.clades.append(clade)
            del clade.parent
            return parent
       
