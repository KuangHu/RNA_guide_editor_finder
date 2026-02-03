# Header Formats and Space Handling

## Summary Table

| Format | Delimiter | Spaces in IDs? | Notes |
|--------|-----------|----------------|-------|
| Diamond BLASTP | TAB | ✓ YES | Subject IDs can have spaces |
| HMMER tblout | SPACE | ✗ NO | Target/query names cannot have spaces |
| Prodigal FAA | FASTA | ✓ YES | Only first word used as ID |
| GFF/GFF3 | TAB | ✓ YES | Tab-delimited |
| BED | TAB | ✓ YES | Tab-delimited |
| FASTA | FASTA | ✓ YES | Only first word used as ID |

## Detailed Format Specifications

### 1. Diamond BLASTP (.m8 format)

**Delimiter:** TAB (`\t`)
**Spaces in names:** ✓ YES (safe)

**Example:**
```
contig_1_1	IS110 family transposase	95.5	200	9	0	1	200	15	214	1e-100	350.5
```

**Parsing method:**
```python
fields = line.split('\t')  # Tab-delimited
query_id = fields[0]      # Can have spaces (rare)
subject_id = fields[1]    # CAN have spaces (common!)
pident = float(fields[2])
```

**Real-world example:**
```
contig_1_1	IS110 family transposase	95.5	...
contig_1_2	transposase [Escherichia coli]	88.2	...
```

✓ **Our parser handles this correctly** - uses `split('\t')`

---

### 2. HMMER tblout format

**Delimiter:** SPACE (multiple spaces)
**Spaces in names:** ✗ NO (will break parsing!)

**Header:**
```
# target name        accession   tlen query name           accession   qlen   E-value  score  bias
```

**Example:**
```
contig_1_1           -            200 IS110_transposase    PF03400.15   250  1.5e-100  350.5   0.1
```

**Parsing method:**
```python
fields = line.split()     # Splits on ALL whitespace
target_name = fields[0]   # MUST NOT contain spaces
query_name = fields[3]    # MUST NOT contain spaces
evalue = float(fields[6])
```

**⚠️ IMPORTANT LIMITATIONS:**

1. **Target names (protein IDs)** CANNOT contain spaces
   - ✓ Good: `contig_1_1`, `protein_A`, `seq123`
   - ✗ Bad: `contig 1 1`, `protein A`, `seq 123`

2. **Query names (HMM profile names)** CANNOT contain spaces
   - ✓ Good: `IS110_transposase`, `PF03400`
   - ✗ Bad: `IS110 transposase`, `PF 03400`

3. **Why?** HMMER uses space-delimited format with fixed-width columns

**What breaks:**
```
# This line will fail to parse correctly:
contig 1 1           -            200 IS110 transposase    PF03400.15 ...

# Fields will be:
# [0] = 'contig'      (should be 'contig 1 1')
# [1] = '1'           (should be '-')
# [2] = '1'           (should be '200')
# ... everything shifted!
```

⚠️ **Our parser:** Uses `split()` - assumes no spaces in names (standard HMMER format)

---

### 3. Prodigal FAA Headers

**Format:** FASTA with special comment structure
**Spaces in protein_id:** ✓ YES (but only first word used)

**Example:**
```
>contig_1_1 # 100 # 300 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=None
MKLVPQRSTAVILGKLMNPQ
```

**Parsing method:**
```python
protein_id = line[1:].split()[0]  # Everything after '>' until first space
parts = line.split('#')
start = int(parts[1].strip())
end = int(parts[2].strip())
strand = int(parts[3].strip())
```

**What gets extracted:**
```
Header: >contig_1_1 # 100 # 300 # 1 # ID=1_1...
protein_id: 'contig_1_1'  (everything before first space)
```

✓ **Our parser handles this correctly** - only uses first field

---

### 4. GFF/GFF3 Files

**Delimiter:** TAB (`\t`)
**Spaces in fields:** ✓ YES (safe)

**Example:**
```
contig_1	Prodigal	CDS	100	300	.	+	0	ID=gene1;product=hypothetical protein
```

**Parsing method:**
```python
fields = line.split('\t')  # Tab-delimited
seqid = fields[0]          # Can have spaces
feature_type = fields[2]   # Can have spaces
attributes = fields[8]     # CAN have spaces (e.g., "product=IS110 transposase")
```

✓ **Our parser handles this correctly** - uses `split('\t')`

---

### 5. BED Files

**Delimiter:** TAB (`\t`)
**Spaces in names:** ✓ YES (safe)

**Example:**
```
contig_1	99	300	gene1	100	+
```

**Parsing method:**
```python
fields = line.split('\t')  # Tab-delimited
chrom = fields[0]          # Can have spaces
name = fields[3]           # Can have spaces
```

✓ **Our parser handles this correctly** - uses `split('\t')`

---

### 6. FASTA Files

**Format:** Standard FASTA
**Spaces in ID:** ✓ YES (but only first word used as ID)

**Example:**
```
>contig_1 Additional description text here
ATCGATCGATCGATCG
```

**Parsing method:**
```python
seq_id = line[1:].strip().split()[0]  # First word after '>'
```

**What gets extracted:**
```
Header: >contig_1 Additional description
ID: 'contig_1'
Description: 'Additional description'
```

✓ **Our parser handles this correctly** - provides both `id` and `description`

---

## Recommendations

### For Your Data

1. **Protein IDs (Prodigal output):**
   - Use underscores instead of spaces: `contig_1_1` ✓
   - Avoid: `contig 1 1` ✗

2. **Contig IDs:**
   - Use underscores: `my_contig_1` ✓
   - Avoid spaces if possible, but our parsers can handle them for most formats

3. **HMM Profile Names:**
   - No spaces allowed: `IS110_transposase` ✓
   - Avoid: `IS110 transposase` ✗

4. **Diamond Subject IDs:**
   - Spaces are OK: `IS110 family transposase` ✓
   - Our parser will preserve them

### Testing Your Data

Use this script to check your headers:

```python
# Check if your HMMER output has spaces in names
with open("hmmer_output.tbl") as f:
    for line in f:
        if not line.startswith('#') and line.strip():
            fields = line.split()
            target_name = fields[0]
            if ' ' in target_name:
                print(f"WARNING: Space in target name: {target_name}")

# Check protein IDs from Prodigal
with open("proteins.faa") as f:
    for line in f:
        if line.startswith('>'):
            protein_id = line[1:].split()[0]
            full_id = line[1:].split('#')[0].strip()
            if protein_id != full_id:
                print(f"Note: ID shortened from '{full_id}' to '{protein_id}'")
```

---

## Quick Reference

**When spaces are safe:**
- ✓ Diamond BLASTP subject IDs
- ✓ GFF attributes
- ✓ BED names
- ✓ FASTA descriptions (after ID)

**When spaces break parsing:**
- ✗ HMMER target names (protein IDs)
- ✗ HMMER query names (HMM profiles)

**When spaces are ignored:**
- ⚠️ Prodigal protein IDs (only first word used)
- ⚠️ FASTA IDs (only first word used)
