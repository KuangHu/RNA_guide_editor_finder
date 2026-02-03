# Interactive Parser Tutorial - Jupyter Notebook Guide

## Quick Start

### Option 1: Using the launch script
```bash
# From the project directory
./launch_notebook.sh
```

### Option 2: Manual launch
```bash
# Activate conda environment
conda activate opfi

# Install jupyter (if not already installed)
pip install jupyter notebook

# Launch notebook
jupyter notebook Interactive_Parser_Tutorial.ipynb
```

### Option 3: Using JupyterLab (recommended)
```bash
conda activate opfi
pip install jupyterlab
jupyter lab Interactive_Parser_Tutorial.ipynb
```

## What's in the Notebook?

The Interactive Parser Tutorial notebook contains **9 comprehensive sections**:

### Part 1: Setup
- Import all modules
- Create sample data files automatically
- No external files needed!

### Part 2: Prodigal Parsers
- `parse_prodigal_faa()` - Parse positions only
- `parse_prodigal_faa_full()` - Parse with sequences
- See input files and outputs interactively

### Part 3: FASTA Parsers
- `parse_fasta()` - Parse as dictionary
- `parse_fasta_records()` - Parse as record objects
- View actual sequences

### Part 4: Diamond BLASTP Parsers
- `parse_diamond_blastp()` - Basic parsing
- `parse_diamond_blastp_with_positions()` - With genomic positions
- See all DiamondBlastHit fields

### Part 5: HMMER Parsers
- `parse_hmm_tblout()` - Basic parsing
- `parse_hmm_tblout_with_positions()` - With genomic positions
- See all HmmHit fields

### Part 6: Sequence Retrieval
- `reverse_complement()` - Reverse complement DNA
- `retrieve_sequence_from_fasta()` - Get specific regions
- `retrieve_sequences_from_hits()` - Batch retrieval
- See sequences added to hits

### Part 7: DataFrame Conversion
- `diamond_hits_to_dataframe()` - Convert to pandas
- `hmm_hits_to_dataframe()` - Convert to pandas
- Explore DataFrame structure

### Part 8: DataFrame Analysis
- Filter by E-value, identity, strand
- Sort by different columns
- Group by subject
- Calculate statistics
- Value counts and aggregations

### Part 9: Export Results
- `save_hits_to_csv()` - Export to CSV
- `df.to_csv()` - Direct pandas export
- Read back and verify

### Part 10: Complete Workflow
- End-to-end example
- All steps combined
- See the complete pipeline

## Features

✅ **No external files needed** - Creates sample data automatically
✅ **Interactive** - Run each cell to see outputs
✅ **Explanatory** - Shows inputs and outputs clearly
✅ **Comprehensive** - Covers all major functions
✅ **Practical** - Real-world examples
✅ **Educational** - Learn by doing

## How to Use

1. **Read through sequentially** - Start from the top
2. **Run each cell** - Press Shift+Enter or click "Run"
3. **Observe outputs** - See what each function returns
4. **Experiment** - Modify parameters and re-run
5. **Try your own data** - Replace sample files with your data

## Tips

### Cell Execution
- **Run a cell**: Shift+Enter
- **Run and stay**: Ctrl+Enter
- **Insert cell below**: B (in command mode)
- **Delete cell**: DD (in command mode)

### Modify Examples
```python
# Change parameters
position_map = parse_prodigal_faa("YOUR_FILE.faa")

# Try different filters
high_conf = df[(df['evalue'] < 1e-100) & (df['pident'] > 95)]

# Add your own analysis
df.groupby('contig_id')['pident'].mean()
```

### Visualization
Add this cell to create plots:
```python
import matplotlib.pyplot as plt

# Plot identity distribution
df['pident'].hist(bins=20)
plt.xlabel('Percent Identity')
plt.ylabel('Count')
plt.title('Identity Distribution')
plt.show()

# Plot E-value vs Identity
plt.scatter(df['evalue'], df['pident'])
plt.xscale('log')
plt.xlabel('E-value (log scale)')
plt.ylabel('Percent Identity')
plt.show()
```

## Troubleshooting

### Jupyter not installed
```bash
conda activate opfi
pip install jupyter notebook
```

### Cannot find module
```bash
# Make sure you're in the project directory
cd /home/kuangh/scripts/IS110/tools/RNA_guide_editor_finder

# Check path
import sys
print(sys.path)
```

### Kernel dies
- Try restarting kernel: Kernel → Restart
- Check memory usage
- Try running cells individually

### Import errors
```bash
# Reinstall dependencies
conda activate opfi
pip install pandas
```

## Example Output

When you run the notebook, you'll see outputs like:

```
INPUT FILE (proteins.faa):
======================================================================
>contig_1_1 # 100 # 300 # 1 # ID=1_1;partial=00;start_type=ATG...
MKLVPQRSTAVILGKLMNPQRSTAVILGKLMNPQ
...
======================================================================

FUNCTION: parse_prodigal_faa()
======================================================================

OUTPUT (position_map):
Type: <class 'dict'>
Number of proteins: 3

Contents:
  contig_1_1: 100..300 (forward)
  contig_1_2: 500..800 (reverse)
  contig_2_1: 50..250 (forward)
```

## Next Steps After the Notebook

Once you're comfortable with the tutorial:

1. **Use with your own data**
   ```python
   position_map = parse_prodigal_faa("/path/to/your/genome.faa")
   hits = parse_diamond_blastp_with_positions("your_results.m8", position_map)
   ```

2. **Create analysis scripts**
   - Copy cells into .py files
   - Automate workflows
   - Process multiple genomes

3. **Advanced analysis**
   - Combine with other tools
   - Create visualizations
   - Generate reports

## Additional Resources

- **README.md** - Main project documentation
- **DATAFRAME_GUIDE.md** - Detailed DataFrame operations
- **HEADER_FORMATS.md** - File format specifications
- **PARSER_QUICK_REFERENCE.md** - Quick function reference
- **examples/demo_*.py** - Standalone demo scripts

## Getting Help

If you encounter issues:

1. Check the error message
2. Review the relevant documentation
3. Try the demo scripts: `python examples/demo_parsers.py`
4. Run tests: `python -m pytest tests/test_parsers.py -v`

## Saving Your Work

The notebook will save automatically, but you can also:

```bash
# Save as Python script
jupyter nbconvert --to python Interactive_Parser_Tutorial.ipynb

# Save as HTML
jupyter nbconvert --to html Interactive_Parser_Tutorial.ipynb

# Save as PDF (requires LaTeX)
jupyter nbconvert --to pdf Interactive_Parser_Tutorial.ipynb
```

---

**Enjoy exploring the parser modules interactively!** 🎉
