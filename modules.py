from io import StringIO, TextIOBase
import re

def load_fasta(fafile):
    entries = {}
    current_header = None
    current_sequence = []
    fh=None
    if  isinstance(fafile, TextIOBase):
        fh=fafile
    else:
        fh=open(fafile,'r')
    
    for line in fh:
        line = line.rstrip()  # remove trailing \n
        
        if line.startswith('>'):
            # Save previous entry if exists
            if current_header is not None:
                seq = ''.join(current_sequence)
                entries[current_header]=seq
            
            # Start new entry; get rid of >
            current_header = line[1:]
            # split and get rid of anything after whitespace
            if re.search(r'\s', current_header):
                current_header=re.split(r'\s+', current_header)[0]
            current_sequence = []
        elif line and current_header is not None:
            # Add sequence line (skip empty lines)
            current_sequence.append(line)
    
    # Don't forget the last entry!
    if current_header is not None:
        seq = ''.join(current_sequence)
        entries[current_header]=seq
    fh.close()
    
    return entries








def test_fasta_loader():
    from io import StringIO
    
    test_content = """>seq1 some description
ACGTACGT
GCTA
>seq2
NNNNNNNNNN
>seq3 empty sequence

>seq4
ATGCATGCATGC
"""

    result = load_fasta(StringIO(test_content))
    
    expected = {
        "seq1": "ACGTACGTGCTA",
        "seq2": "NNNNNNNNNN",
        "seq3": "",
        "seq4": "ATGCATGCATGC"
    }
    
    assert len(result) == 4,           f"Expected 4 sequences, got {len(result)}"
    assert "seq3" in result,           "Missing empty sequence entry"
    assert result["seq3"] == "",       "Empty sequence should be empty string"
    assert result == expected,         "Dictionary content doesn't match expected"
    
    print("Quick test of fasta_loader PASSED ✓")

if __name__ == "__main__":
    test_fasta_loader()