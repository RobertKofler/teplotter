#!/bin/bash

set -euo pipefail

mkdir -p "$PREFIX/bin"
mkdir -p "$PREFIX/lib/seqvista"

cp src/*.py "$PREFIX/lib/seqvista/"
cp src/*.R  "$PREFIX/lib/seqvista/"

# Create the SeqVista entry point
cat > "$PREFIX/bin/SeqVista" << 'EOF'
#!/bin/bash
exec python "$CONDA_PREFIX/lib/seqvista/seqvista.py" "$@"
EOF
chmod +x "$PREFIX/bin/SeqVista"
