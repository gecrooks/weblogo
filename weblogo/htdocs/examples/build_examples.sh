#!/bin/bash

cd "$(dirname "${BASH_SOURCE[0]}")"

echo -n '.'
uv run weblogo --format png --size large \
    -i -5 -l 1 -u 20 \
    --title "The DNA-binding helix-turn-helix motif of the CAP family" \
    <   cap_hth.fa    \
    > cap_hth.png || exit
    
echo -n '.'
uv run weblogo --format png --size large \
    --title "58 CAP Binding Sites" \
    -i -10                         \
    < cap_dna.fa    \
    > cap_dna.png || exit


echo -n '.'
uv run weblogo --format png --size large \
    --title "19 LexA Binding Sites" \
    -i -9                           \
    < lexA.fa    \
    > lexA.png    || exit
    
echo -n '.'
uv run weblogo --format png --size large \
    --title "-10 region of E. coli promoters" \
    -i -21 --lower 0 -u 7               \
    < ecoli10.fa    \
    > ecoli10.png || exit
    
echo -n '.'
uv run weblogo --format png --size large \
    -l 63 -u 83  \
    < globins.fa    \
    > globins.png || exit

#echo -n '.'
#weblogo --format PNG --size large \
#    -l 31 -u 150  \
#    < globins.fa    \
#    > more_globins.png || exit
    
echo -n '.'
uv run weblogo --format png --size large \
    --title "Helix-Turn-Helix Motifs" \
    -i -11 -l 1 -u 17  --yaxis 3.5 \
    < hth.fa    \
    > hth.png || exit

echo -n '.'
uv run weblogo --format png --size large \
    --title "exon|intron" \
    -i -11 -l -6 -u 8 \
    < exon-intron.fa    \
    > exon-intron.png || exit

echo -n '.'
uv run weblogo --format png --size large \
    --title "intron | exon" \
    -i -21 -l -20 -u 3 \
    < intron-exon.fa    \
    > intron-exon.png || exit

echo

# --- SVG examples (native renderer) ---
echo -n ' svg:'
echo -n '.'
uv run weblogo --format svg --size large \
    -i -5 -l 1 -u 20 \
    --title "The DNA-binding helix-turn-helix motif of the CAP family" \
    < cap_hth.fa > cap_hth.svg || exit

echo -n '.'
uv run weblogo --format svg --size large \
    --title "58 CAP Binding Sites" \
    -i -10                         \
    < cap_dna.fa > cap_dna.svg || exit

echo -n '.'
uv run weblogo --format svg --size large \
    -l 63 -u 83  \
    < globins.fa > globins.svg || exit

echo

# --- PDF examples (native renderer) ---
echo -n ' pdf:'
echo -n '.'
uv run weblogo --format pdf --size large \
    -i -5 -l 1 -u 20 \
    --title "The DNA-binding helix-turn-helix motif of the CAP family" \
    < cap_hth.fa > cap_hth.pdf || exit

echo -n '.'
uv run weblogo --format pdf --size large \
    --title "58 CAP Binding Sites" \
    -i -10                         \
    < cap_dna.fa > cap_dna.pdf || exit

echo -n '.'
uv run weblogo --format pdf --size large \
    -l 63 -u 83  \
    < globins.fa > globins.pdf || exit

echo

