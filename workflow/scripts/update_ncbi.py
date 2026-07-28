#!/usr/bin/env python3
import sys
from ete3 import NCBITaxa
from pathlib import Path

db_dir = Path(sys.argv[1])
db_dir.mkdir(parents=True, exist_ok=True)
db_file = db_dir / "taxa.sqlite"

ncbi = NCBITaxa(dbfile=str(db_file))
ncbi.update_taxonomy_database()