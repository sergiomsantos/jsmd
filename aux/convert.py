import xmltodict

def convert(xml_file, xml_attribs=True):
    with open(xml_file) as f:
        d = xmltodict.parse(f, xml_attribs=xml_attribs)
        return d

if __name__ == '__main__':
    import json
    from sys import argv

    data = convert(argv[1])
    print json.dumps(data, indent=2).replace('@','')