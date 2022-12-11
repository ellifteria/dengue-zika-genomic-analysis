# gs -dSAFER -dQUIET -dNOPLATFONTS -dNOPAUSE -dBATCH \ 
#   # When converting multiple-page PDFs you should add "%d" to the filename-string 
#   # which will be replaced with the sequence-number of the file
#   -sOutputFile="$2" \ 
#   # Create a PNG-File with alpha-transparency
#   -sDEVICE=pngalpha
#   # resolution in DPI
#   -r72 \ 
#   # Use high grade text antialiasing. Can be 0, 1, 2 or 4
#   -dTextAlphaBits=4 \
#   # Use high grade graphics anti-aliasing. Can be 0, 1, 2 or 4
#   -dGraphicsAlphaBits=4 \
#   # If you are converting a CMYK-PFD to RGB-color you should use CIE-Color
#   -dUseCIEColor \
#   # use the PDFs Trim-Box to define the final image
#   -dUseTrimBox \
#   # start converting on the first side of the PDF
#   -dFirstPage=1 \
#   # convert only until the first page of the PDF
#   -dLastPage=1 \
#   # Path to the file you want to convert
#   $1

convert -density 600 $1 -quality 100 -fill white -opaque none -alpha off $2_%d.png