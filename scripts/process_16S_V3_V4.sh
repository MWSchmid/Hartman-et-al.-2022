#!/bin/bash
#
# process_16S_V3_V4.sh
#
# Authors: Marc W. Schmid <contact@mwschmid.ch>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

me=$(basename "$0")
curPath=$(dirname "$0")

## defaults
rawReadsDir=""
processedReadsDir=""
usearchOutput=""
utaxDir=""
scriptDir=""
tempDir=""
threads=1
memory=4
maxMem=4

PRIMER_16S_FORWARD="CCTACGGGNGGCWGCAG"
PRIMER_16S_REVERSE="GACTACHVGGGTATCTAATCC"

## Exit status codes (following <sysexits.h>)
EX_OK=0			# successful termination
EX_USAGE=64		# command line usage error
EX_DATAERR=65		# data format error
EX_NOINPUT=66		# cannot open input
EX_NOUSER=67		# addressee unknown
EX_NOHOST=68		# host name unknown
EX_UNAVAILABLE=69	# service unavailable
EX_SOFTWARE=70		# internal software error
EX_OSERR=71		# system error (e.g., can't fork)
EX_OSFILE=72		# critical OS file missing
EX_CANTCREAT=73		# can't create (user) output file
EX_IOERR=74		# input/output error
EX_TEMPFAIL=75		# temp failure; user is invited to retry
EX_PROTOCOL=76		# remote error in protocol
EX_NOPERM=77		# permission denied
EX_CONFIG=78		# configuration error

## helper functions

function die () {
    rc="$1"
    shift
    (echo -n "$me: ERROR: ";
        if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
    exit $rc
}

## usage info

usage () {
    cat <<__EOF__
Usage:
  $me [options] rawReadsDir processedReadsDir usearchOutput utaxDir scriptDir tempDir
  
The script will run:
1) Fastp to remove adapters
2) usearch -fastq_mergepairs to merge the reads
3) trim16Smanual.py to remove everything outside of the 16S primers
4) usearch -fastq_filter to filter for high quality reads (requires file splitting)
5) fqtrim to dereplicate sequences (instead of usearch due to 32 bit limitations)
6) extract sequences that were seen at least twice and sort
7) usearch -unoise3 do get zero-OTUs
8) usearch -cluster_smallmem to get the 99% OTUs
9) usearch -sintax to predict the taxonomy
10) usearch -otutab to count reads per OTU (requires file splitting)
11) file compression
  
NOTE:
The primer sequences are hard-coded:

forward: CCTACGGGNGGCWGCAG
reverse: GACTACHVGGGTATCTAATCC
  
Arguments:
rawReadsDir: Directory with the original read files. Fastq files must be named SampleX_R1/2.fastq.gz.
processedReadsDir Directory in which all processed read files will be stored.
usearchOutput: Directory for the Usearch output.
utaxDir: Directory with the RDP database is located (rdp_16s_v16_sp.fa).
scriptDir: Directory with all required scripts (trim16Smanual.py, pUPARSE_extractMinRepCountSeqs.py, pUPARSE_sortDerep.py, pUPARSE_splitOneFastaIntoSeveral.py, mergeCountTables.R).
tempDir: A temporary directory.

Options:
  -v            Enable verbose logging (no effect)
  -h            Print this help text
  -t		    Number of available threads
  -m            Amount of memory to be allocated (per core, in GB)
__EOF__
}

warn () {
  (echo -n "$me: WARNING: ";
      if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
}

have_command () {
  type "$1" >/dev/null 2>/dev/null
}

require_command () {
  if ! have_command "$1"; then
    die $EX_UNAVAILABLE "Could not find required command '$1' in system PATH. Aborting."
  fi
}

is_absolute_path () {
    expr match "$1" '/' >/dev/null 2>/dev/null
}

input_exists () {
  echo -n "Checking input file ${1}..."
  if [ -e $1 ]; then
    echo "[ok]"
  else
    echo "[failed]"
    die $EX_NOINPUT "Could not find input file: ${1}."
  fi
}

output_exists () {
  echo -n "Checking output file ${1}... "
  if [ -e $1 ]; then
    echo "[ok]"
  else
    echo "[failed]"
    die $EX_OSFILE "Could not find output file: ${1}."
  fi
}

remove_if_present () {
  if [ -e $1 ]; then
    rm $1
  fi
}

remove_dir_if_present () {
  if [ -e $1 ]; then
    rm -r $1
  fi
}

## parse command-line

short_opts='hvt:m:'
long_opts='help,verbose,threads,memory'

getopt -T > /dev/null
rc=$?
if [ "$rc" -eq 4 ]; then
    # GNU getopt
    args=$(getopt --name "$me" --shell sh -l "$long_opts" -o "$short_opts" -- "$@")
    if [ $? -ne 0 ]; then
        die $EX_USAGE "Type '$me --help' to get usage information."
    fi
    # use 'eval' to remove getopt quoting
    eval set -- $args
else
    # old-style getopt, use compatibility syntax
    args=$(getopt "$short_opts" "$@")
    if [ $? -ne 0 ]; then
        die $EX_USAGE "Type '$me --help' to get usage information."
    fi
    set -- $args
fi

while [ $# -gt 0 ]; do
    case "$1" in
	--threads|-t)  		shift; threads=$1 ;;
	--memory|-m)  		shift; memory=$1 ;;	
    --verbose|-v) 		verbose='--verbose' ;;
    --help|-h)    		usage; exit 0 ;;
    --)           		shift; break ;;
    esac
    shift
done

maxMem=$(($threads * $memory))

## sanity checks

# the required arguments must be present 
if [ $# -lt 5 ]; then
    die $EX_USAGE "Missing required arguments. Type '$me --help' to get usage help."
fi

rawReadsDir=$1
shift
processedReadsDir=$1
shift
usearchOutput=$1
shift
utaxDir=$1
shift
scriptDir=$1
shift
tempDir=$1
shift

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

# checking programs and scripts

require_command python2.7
require_command usearch
require_command fastp
require_command fqtrim
require_command Rscript
input_exists "$scriptDir/trim16Smanual.py"
input_exists "$utaxDir/rdp_16s_v16_sp.fa"
input_exists "$scriptDir/pUPARSE_extractMinRepCountSeqs.py"
input_exists "$scriptDir/pUPARSE_sortDerep.py"
input_exists "$scriptDir/pUPARSE_splitOneFastaIntoSeveral.py"
input_exists "$scriptDir/mergeCountTables.R"

###############################################################################################
## Processing begin
###############################################################################################

# Remove Illumina adapters and do some quality trimming with fastp
for curFile in $rawReadsDir/*_R1.fastq.gz; do
inputFileForward=$(basename $curFile)
inputFileReverse=${inputFileForward/_R1/_R2}
prefix=${inputFileForward//_R1.fastq.gz}
command="fastp -w ${threads} -l 30 -i ${rawReadsDir}/${inputFileForward} -I ${rawReadsDir}/${inputFileReverse} -o ${processedReadsDir}/${inputFileForward/.fastq.gz/.tr.fq.gz} -O ${processedReadsDir}/${inputFileReverse/.fastq.gz/.tr.fq.gz} -j ${tempDir}/$prefix.json -h ${tempDir}/$prefix.html"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
mv ${tempDir}/$prefix.json $processedReadsDir/
mv ${tempDir}/$prefix.html $processedReadsDir/
done

# Merge reads
rm $tempDir/*.tr.fq
mv $processedReadsDir/*.tr.fq.gz $tempDir/
gunzip $tempDir/*.tr.fq.gz
command="usearch -fastq_mergepairs $tempDir/*_R1.tr.fq -fastqout $tempDir/merged.fq -relabel @ -report $usearchOutput/report.txt -fastq_maxdiffs 25"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $tempDir/merged.fq
rm $tempDir/*.tr.fq

# Trim manually everything up to the 16S primers, also remove the primers
command="python $scriptDir/trim16Smanual.py 2 10 true $PRIMER_16S_FORWARD $PRIMER_16S_REVERSE $tempDir/merged.fq > $usearchOutput/merged.tr.fq"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $usearchOutput/merged.tr.fq
remove_if_present $tempDir/merged.fq

# Usearch 32 bit has some limitations, we need to split, filter, merge
awk -v TEMPDIR="$tempDir" '{print $0 > TEMPDIR"/tempSeparate_"int((NR-1)/10000000)".fq"}' $usearchOutput/merged.tr.fq
for curFile in $tempDir/tempSeparate_*; do
command="usearch -fastq_filter $curFile -fastaout $curFile.filt -fastq_maxee 1"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
done
cat $tempDir/*.filt > $usearchOutput/filtered.fa
output_exists $usearchOutput/filtered.fa
rm $tempDir/*.filt
rm $tempDir/tempSeparate_*

# Usearch 32 bit has some limitations, we need to use fqtrim for dereplication
command="fqtrim -C -r $usearchOutput/fqtrimReport.txt $usearchOutput/filtered.fa | sed 's/_x/;size=/' > $usearchOutput/derep.fa"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $usearchOutput/derep.fa
output_exists $usearchOutput/fqtrimReport.txt

# Keep only sequences that were seen at least twice
command="python $scriptDir/pUPARSE_extractMinRepCountSeqs.py 2 $usearchOutput/derep.fa $usearchOutput/derep_minRep2.fa"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $usearchOutput/derep_minRep2.fa

# Sort selected and derep sequences (required for next steps)
command="python $scriptDir/pUPARSE_sortDerep.py $usearchOutput/derep_minRep2.fa > $usearchOutput/derep_minRep2_sorted.fa"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $usearchOutput/derep_minRep2_sorted.fa

# Denoise sequences (get z-otus), default minsize of 8
command="usearch -unoise3 $usearchOutput/derep_minRep2_sorted.fa -minsize 8 -zotus $tempDir/zotus.fa"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
sed "s/Zotu/zot/g" $tempDir/zotus.fa > $usearchOutput/zotus.fa
output_exists $usearchOutput/zotus.fa
remove_if_present $tempDir/zotus.fa

# Denoise sequences (get z-otus), reduced minsize of 4
command="usearch -unoise3 $usearchOutput/derep_minRep2_sorted.fa -minsize 4 -zotus $tempDir/zotus_ms4.fa"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
sed "s/Zotu/zot/g" $tempDir/zotus_ms4.fa > $usearchOutput/zotus_ms4.fa
output_exists $usearchOutput/zotus_ms4.fa
remove_if_present $tempDir/zotus_ms4.fa

# Cluster zero-OTUs to get 99 percent OTUs
command="usearch -sortbylength $usearchOutput/zotus.fa -fastaout $usearchOutput/zotus_sorted.fa -minseqlength 64"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
command="usearch -cluster_smallmem $usearchOutput/zotus_sorted.fa -id 0.99 -centroids $usearchOutput/otus.fa"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $usearchOutput/otus.fa

command="usearch -sortbylength $usearchOutput/zotus_ms4.fa -fastaout $usearchOutput/zotus_ms4_sorted.fa -minseqlength 64"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
command="usearch -cluster_smallmem $usearchOutput/zotus_ms4_sorted.fa -id 0.99 -centroids $usearchOutput/otus_ms4.fa"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $usearchOutput/otus_ms4.fa

# Predict taxonomy
command="usearch -sintax $usearchOutput/otus.fa -db $utaxDir/rdp_16s_v16_sp.fa -strand both -tabbedout $usearchOutput/otus_tax.txt -sintax_cutoff 0.8"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
command="usearch -sintax $usearchOutput/zotus.fa -db $utaxDir/rdp_16s_v16_sp.fa -strand both -tabbedout $usearchOutput/zotus_tax.txt -sintax_cutoff 0.8"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $usearchOutput/otus_tax.txt
output_exists $usearchOutput/zotus_tax.txt

command="usearch -sintax $usearchOutput/otus_ms4.fa -db $utaxDir/rdp_16s_v16_sp.fa -strand both -tabbedout $usearchOutput/otus_ms4_tax.txt -sintax_cutoff 0.8"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
command="usearch -sintax $usearchOutput/zotus_ms4.fa -db $utaxDir/rdp_16s_v16_sp.fa -strand both -tabbedout $usearchOutput/zotus_ms4_tax.txt -sintax_cutoff 0.8"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
output_exists $usearchOutput/otus_ms4_tax.txt
output_exists $usearchOutput/zotus_ms4_tax.txt

# Count reads, again split and merge
command="python $scriptDir/pUPARSE_splitOneFastaIntoSeveral.py $usearchOutput/filtered.fa $tempDir 5000000"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
for curFile in $tempDir/temp_*.fasta; do
command="usearch -otutab $curFile -otus $usearchOutput/otus.fa -strand plus -otutabout ${curFile//.fasta}.otu.counts"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
command="usearch -otutab $curFile -zotus $usearchOutput/zotus.fa -strand plus -otutabout ${curFile//.fasta}.zotu.counts"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
command="usearch -otutab $curFile -otus $usearchOutput/otus_ms4.fa -strand plus -otutabout ${curFile//.fasta}.otu_ms4.counts"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
command="usearch -otutab $curFile -zotus $usearchOutput/zotus_ms4.fa -strand plus -otutabout ${curFile//.fasta}.zotu_ms4.counts"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
done

# Merge the count tables with an Rscript
command="Rscript $scriptDir/mergeCountTables.R $tempDir $usearchOutput"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"

# Delete temp files
rm $tempDir/temp_*

# Compress some files
for curFile in $usearchOutput/*.fa $usearchOutput/*.fq; do
command="pigz -p $threads $curFile"
echo "=== ${me}: Running: ${command}"; eval $command; rc=$?; echo "=== ${me}: Command ended with exit code $rc"
done

###############################################################################################
## Processing end
###############################################################################################

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."
exit $rc





