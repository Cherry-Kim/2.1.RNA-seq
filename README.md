# RNA-seq
# Download the latest release from and uncompress it
$wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz

$tar -xzf 2.7.9a.tar.gz

# Compile under Linux
$cd STAR-2.7.9a/source/

$make STAR

$sudo cp STAR /usr/local/bin
