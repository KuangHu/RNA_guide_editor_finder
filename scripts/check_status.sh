#!/bin/bash
# Status checker for GTDB indexing and search

echo "========================================"
echo "GTDB Processing Status Check"
echo "========================================"
echo ""

# Check indexing process
echo "[1] Indexing Process:"
if ps aux | grep -q "mmseqs createindex" | grep -v grep; then
    echo "  ✓ RUNNING"
    ps aux | grep "mmseqs createindex" | grep -v grep | awk '{print "    PID:", $2, "CPU:", $3"%", "MEM:", $4"%"}'
else
    echo "  ✗ Not running"
fi
echo ""

# Check index directory size
echo "[2] Index Directory Size:"
if [ -d "/groups/rubin/projects/kuang/db/GTDB_organized/mmseqs_db/tmp_index" ]; then
    du -sh /groups/rubin/projects/kuang/db/GTDB_organized/mmseqs_db/tmp_index/
else
    echo "  tmp_index/ does not exist"
fi
echo ""

# Check for .idx file
echo "[3] Index File (.idx):"
if [ -f "/groups/rubin/projects/kuang/db/GTDB_organized/mmseqs_db/gtdb_mmseq_db.idx" ]; then
    ls -lh /groups/rubin/projects/kuang/db/GTDB_organized/mmseqs_db/gtdb_mmseq_db.idx
    echo "  ✓ INDEX READY!"
else
    echo "  Not created yet (still indexing)"
fi
echo ""

# Check log
echo "[4] Latest Index Log (last 10 lines):"
if [ -f "/tmp/claude/-home-kuangh-scripts-IS110-tools-RNA-guide-editor-finder/tasks/bc4d301.output" ]; then
    tail -10 /tmp/claude/-home-kuangh-scripts-IS110-tools-RNA-guide-editor-finder/tasks/bc4d301.output
else
    echo "  No log found"
fi
echo ""

# System resources
echo "[5] System Resources:"
echo "  CPU cores: $(nproc)"
free -h | grep -E "Mem:|Swap:"
echo ""

# Directory organization
echo "[6] GTDB Organization:"
ls -lh /groups/rubin/projects/kuang/db/GTDB_organized/
echo ""

echo "========================================"
echo "To monitor in real-time:"
echo "  watch -n 5 $0"
echo "========================================"
