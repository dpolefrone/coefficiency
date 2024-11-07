from pypdf import PdfReader

# initial code AFZ 110424; DAP/AMA ongoing edits

##!!!! I didn't realize these were all in ethanol, you may need to update this list for different solvents!!!##

common_solvents = [
    "ethanol", "methanol", "water", "acetone", "acetonitrile", "dichloromethane",
    "chloroform", "toluene", "hexane", "cyclohexane", "isopropanol", "THF", "DMF", "DMSO"]

# set of 38 items (378.8MB) in ./reaxysRN_records; block[A..AL].pdf. 
pdf_path = "./reaxysRN_records/blockA.pdf"  # Replace with the path to your PDF file
pdf_reader = PdfReader(pdf_path)

text = ""
for page in pdf_reader.pages:
    text += page.extract_text()

lines = text.splitlines()
rows = []
current_row = []
for line in lines:
    if "Reaxys ID" in line and current_row:
        rows.append(current_row)
        current_row = []
    current_row.append(line)

if current_row:
    rows.append(current_row)

def extract_details(row):
    cas_number = None
    inchi_key = None
    entries = []

    for i, line in enumerate(row):
        print(line)

        '''
        if "CAS Registry Number" in line:
            cas_number = line.split(":")[-1].strip()
        if "InChi Key" in line: 
            inchi_key = line.split(":")[-1].strip()

        if "Ext./Abs. Coefficient" in line:
            ext_coeff = row[i + 1].strip() if i + 1 < len(row) else None
            abs_maxima = None
            solvent = None
            for j in range(i - 3, -1, -1):  
                if "Absorption Maxima" in row[j].lower(): 
                    abs_maxima = row[i + 1].strip() if i + 1 < len(row) else None
                if any(solvent_name.lower() in row[j].lower() for solvent_name in common_solvents):
                    solvent = row[j].strip()
                else: 
                    solvent = None
            if ext_coeff:
                entries.append((ext_coeff, abs_maxima, solvent))
                print(entries)

    return cas_number, entries
    '''

extracted_data = []
for row in rows:
    cas_number, entries = extract_details(row)
    
    for ext_coeff, abs_maxima, solvent in entries:
        if cas_number and ext_coeff and abs_maxima and solvent:
            extracted_data.append({
                "CAS Number": cas_number,
                "Extinction Coefficient": ext_coeff.strip('Copyright © 2024 Elsevier Life Sciences IP Limited except certain content provided by'),
                "Absorbance Maxima": abs_maxima,
                "Solvent": solvent.strip('scopy)')
            })

for data in extracted_data:
    print(data)
