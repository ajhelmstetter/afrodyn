#!/bin/bash

# After alignment this script must be run in the output folder
# moves trimmed alignments to their own folder "gblocks" and renames them
# for downstream use

mkdir gblocks

mv *-gb gblocks

cd gblocks

rename -v 'A-gb' 'A' *

