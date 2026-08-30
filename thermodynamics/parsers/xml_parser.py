import xml.etree.ElementTree as ET


class XMLParser:

    def __init__(self, xml_file):
        self.xml_file = xml_file
        self.tree = ET.parse(xml_file)
        self.root = self.tree.getroot()


    def _local_name(self, tag):
        """Extract the local name from an XML tag."""
        if '}' in tag:
            return tag.split('}', 1)[1]
        return tag

    def _element_text(self, element):
        """Extract the text content from an XML element."""
        return element.text.strip() if element.text else ''

    def _normalize_text(self, text):
        """Normalize text by removing extra whitespace and newlines."""
        return ' '.join(text.split())