#!/bin/bash

mkdir -p tmp

echo "# Test weblogo by building logos with many different options."

echo  -n '.'
./weblogo  < cap.fa > tmp/logo0.eps ||exit

echo  -n '.'
./weblogo --title "Default Logo with Title" < cap.fa > tmp/logo1.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint "Default Logo with this fineprint and debug on" < cap.fa > tmp/logo2.eps ||exit

echo  -n '.'
./weblogo --debug no  --fineprint "--debug no"   --debug no < cap.fa > tmp/logo3.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint "" --title "No fine print" --debug yes < cap.fa > tmp/logo4.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint "No title" --title "" < cap.fa > tmp/logo5.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint "--first-index -10" --first-index -10 < cap.fa > tmp/logo6.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint " --first-index -10 --stacks-per-line 11 " --first-index -10 --stacks-per-line 11 < cap.fa > tmp/logo7a.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint " --first-index -10 --stacks-per-line 8 " --first-index -10 --stacks-per-line 8 < cap.fa > tmp/logo7b.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint " --first-index -10 --stacks-per-line 7 " --first-index -10 --stacks-per-line 7 < cap.fa > tmp/logo7c.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint "Test fin and fout" --fin cap.fa  --fout logo8.eps ||exit

# Test Y Axis

echo  -n '.'
./weblogo --debug yes --fineprint "Custom yaxis label " --ylabel 'yaxis label' < cap.fa > tmp/logo9a.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint "Custom units" --units 'nats' < cap.fa > tmp/logo9b.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint "Override custom units with custom yaxis label."  --ylabel 'yaxis label' --units nats < cap.fa > tmp/logo9c.eps ||exit

echo  -n '.'
./weblogo --debug yes --fineprint "Empty ylabel"  --ylabel '' < cap.fa > tmp/logo9d.eps

echo  -n '.'
./weblogo --debug yes --fineprint "No Yaxis"  --show-yaxis no  < cap.fa > tmp/logo9e.eps ||exit

# Test X Axis

echo  -n '.'
./weblogo --debug yes --format pdf --fineprint "Custom xaxis label " --xlabel 'xaxis label' < cap.fa > tmp/logo10a.pdf ||exit

echo  -n '.'
./weblogo --debug yes --format pdf --fineprint "Empty xlabel"  --xlabel '' < cap.fa > tmp/logo10b.pdf ||exit

echo  -n '.'
./weblogo --debug yes --format pdf --fineprint "No Xaxis"  --show-xaxis no  < cap.fa > tmp/logo10c.pdf ||exit

echo  -n '.'
./weblogo --debug yes --format pdf --fineprint "No Xaxis, custom label"  --xlabel "Custom xlabel" --show-xaxis no  < cap.fa > tmp/logo10d.pdf ||exit

# Test Formats

echo  -n '.'
./weblogo --debug no  --fineprint "Format: eps" --format eps < cap.fa > tmp/logo11a.eps ||exit

echo  -n '.'
./weblogo --debug no  --fineprint "Format: png" --size large --format png < cap.fa > tmp/logo11b.png ||exit

echo  -n '.'
./weblogo --debug no  --fineprint "Format: png high res" --format png_print < cap.fa > tmp/logo11c.png ||exit

echo  -n '.'
./weblogo --debug no  --fineprint "Format: pdf" --format pdf < cap.fa > tmp/logo11d.pdf ||exit

echo  -n '.'
./weblogo --debug no  --fineprint "Format: jpeg" --size large --format jpeg < cap.fa > tmp/logo11e.jpeg ||exit

echo  -n '.'
./weblogo --debug no  --fineprint "Format: EPS" --format EPS < cap.fa > tmp/logo11f.eps ||exit

# Test Sizes

echo  -n '.'
./weblogo --debug no  --format png_print --fineprint "default size" < cap.fa > tmp/logo12_default.png ||exit

echo  -n '.'
./weblogo --debug no  --format png_print --fineprint "--size large" --size large < cap.fa > tmp/logo12_large.png ||exit

echo  -n '.'
./weblogo --debug no  --format png_print --fineprint "--size medium" --size medium < cap.fa > tmp/logo12_medium.png ||exit

echo  -n '.'
./weblogo --debug no  --format png_print --fineprint "--size small" --size small < cap.fa > tmp/logo12_small.png ||exit



echo  -n '.'
./weblogo --format pdf --fineprint ""  > tmp/logo13.pdf << LimitString
>
GTTGTTGTTGTT
>
GTCGTCGTCGTC
>
GGGGGGGGGGGG
>
GGAGGAGGAGGA
LimitString




# Test unit options
echo  -n '.'
./weblogo --format pdf --fineprint "probability" --unit probability  > tmp/logo14a.pdf < cap.fa ||exit

echo  -n '.'
./weblogo --format pdf --fineprint "bits" --unit bits  > tmp/logo14b.pdf < cap.fa ||exit

echo  -n '.'
./weblogo --format pdf --fineprint "nats" --unit nats  > tmp/logo14c.pdf < cap.fa ||exit

echo  -n '.'
./weblogo --format pdf --fineprint "kJ/mol" --unit kJ/mol \
     > tmp/logo14d.pdf < cap.fa ||exit

echo  -n '.'
./weblogo --format pdf --fineprint "kT" --unit kT  \
    > tmp/logo14e.pdf < cap.fa ||exit

echo  -n '.'
./weblogo --format pdf --fineprint "kcal/mol" --unit kcal/mol \
    > tmp/logo14f.pdf < cap.fa || exit

echo





