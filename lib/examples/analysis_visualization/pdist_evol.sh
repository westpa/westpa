### SECTION 1 ###

# Run w_pdist first
w_pdist -W west.h5

# Making a folder for the movie
if [ ! -e pdist_evol ];then
  mkdir pdist_evol
else
  # Delete if exists
  rm -rf pdist_evol
  mkdir pdist_evol
  :
fi

### SECTION 2 ###

# We loop over the # of iterations we have and make a plot for each iteration
# This one will be cumulatively done, you can instead make it windowed by 
# changing the --first-iter value around
for i in $(seq 1 1 40);do
  liter=$i
  nfi=$(printf "%04d" $fiter)
  lfi=$(printf "%04d" $liter)
  echo "######## Doing iter $i ########"
  # This will be used by the function in adjust_plot to determine the current 
  # iteration and plot paths
  echo $liter > cur_iter.txt
  plothist average -o pdist_evol/from_${nfi}_${lfi}.png --first-iter 1 --last-iter $liter --postprocess-function postprocess.adjust_plot pdist.h5 0 1
done

### SECTION 3 ###

# Go to the folder with the pngs
cd pdist_evol

# Converting to SGI before encoding
fmt=sgi
for i in *.png
do
  echo $i
  mogrify -format $fmt PNG:$i
done

# Two pass mencoder to encode the movie
mencoder \
  -ovc x264 \
  -x264encopts pass=1:turbo:bitrate=5928:bframes=1:subq=6:frameref=6:me=hex:partitions=all:threads=auto:keyint=300 \
  -mf type=sgi:w=1280:h=720:fps=15 \
  -nosound \
  -o /dev/null mf://\*.sgi

mencoder \
  -ovc x264 \
  -x264encopts pass=2:turbo:bitrate=5928:bframes=1:subq=6:frameref=6:me=hex:partitions=all:threads=auto:keyint=300 \
  -mf type=sgi:w=1280:h=720:fps=15 \
  -nosound \
  -o pdist_evol.avi mf://\*.sgi
rm -f divx2pass.log
