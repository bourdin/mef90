"""Python-Markdown extension turning DOI identifiers into links.

``doi:10.1103/PhysRevLett.87.045501`` becomes a link to
``https://doi.org/10.1103/PhysRevLett.87.045501``.

Used by the FORD documentation build (see ford.md: md_extensions).
"""
import xml.etree.ElementTree as etree

from markdown.extensions import Extension
from markdown.inlinepatterns import InlineProcessor

# 10.NNNN/suffix. The suffix may contain almost anything (including
# parentheses, as in older Elsevier DOIs such as
# 10.1016/S0022-5096(96)00028-X), so stop only at whitespace, quotes,
# and markup characters, and require the final character to be
# alphanumeric so trailing punctuation is not swallowed.
DOI_PATTERN = r"[Dd][Oo][Ii]:\s?(10\.\d{4,9}/[^\s\"'<>\[\]]*[A-Za-z0-9])"


class DoiLinkProcessor(InlineProcessor):
    def handleMatch(self, m, data):
        doi = m.group(1)
        link = etree.Element("a")
        link.set("href", f"https://doi.org/{doi}")
        link.text = f"doi:{doi}"
        return link, m.start(0), m.end(0)


class DoiLinkExtension(Extension):
    def extendMarkdown(self, md):
        # priority between autolink (120) and explicit [text](url)
        # links (160), so hand-written markdown links still win
        md.inlinePatterns.register(DoiLinkProcessor(DOI_PATTERN, md), "doi_link", 140)


def makeExtension(**kwargs):
    return DoiLinkExtension(**kwargs)
