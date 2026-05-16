import re

with open('C:/Users/Ws/WorkBuddy/2026-05-13-task-19/ref_pmc_oai.xml', 'r', encoding='utf-8') as f:
    content = f.read()

# The OAI-PMH format wraps content in XML namespaces
# Find the Full text XML section (the actual JATS XML)
fulltext_start = content.find('<body>')
fulltext_end = content.find('</body>')

if fulltext_start > 0 and fulltext_end > 0:
    body = content[fulltext_start:fulltext_end+7]
    
    # Find Methods section
    methods_start = body.find('<sec')
    methods_match = re.search(r'<title>Methods</title>(.*?)(?=<sec|\Z)', body, re.DOTALL)
    
    if methods_match:
        methods_content = methods_match.group(1)
        # But Methods might be the wrapper and have sub-sections
        # Let's find all content within the Methods sec tag
        text = re.sub(r'<[^>]+>', '\n', methods_content)
        text = text.replace('&#8211;', '-').replace('&amp;', '&').replace('&lt;', '<').replace('&gt;', '>')
        text = re.sub(r'\n\s*\n', '\n\n', text).strip()
        print("=== METHODS CONTENT ===")
        print(text[:10000])
    else:
        # Try to find all text after Abstract
        abstract_end = body.find('</sec>', body.find('</title>', body.find('<title>Abstract')))
        if abstract_end < 0:
            abstract_end = body.find('</sec>', body.find('</title>'))
        rest = body[abstract_end:min(abstract_end+15000, len(body))]
        clean = re.sub(r'<[^>]+>', '\n', rest)
        clean = clean.replace('&#8211;', '-').replace('&amp;', '&')
        clean = re.sub(r'\n\s*\n', '\n\n', clean).strip()
        print("=== AFTER ABSTRACT (first 10000 chars) ===")
        print(clean[:10000])
else:
    print(f"No body found. body at {fulltext_start}")
    # Print part of the file
    print(content[2000:5000])
