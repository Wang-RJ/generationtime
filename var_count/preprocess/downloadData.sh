#!/bin/bash

## This script will create a Data/ directory in the download
## the data used in this study each in its given sub-directory.

shopt -s extglob

# Set the download directory (same as working directory by default):

workingdir="$( pwd )"
outputdir=${workingdir}


# Setting option for downloading only specific dataset (default: all):

declare -a steps2run
steps2run=(step1 step2 step3 step4 step5)
runonlystep=""


# Help menu:

function usage() {
    echo "
        Download Data Help Section:
        ===========================

        Usage: $0
        
        Or for more options:

        Example: $0 -o ~/path/to/my/directory/

        Optional arguments:
        -o|--outputdir <path>       :       Directory path for downloaded data
        --runonlystep <string>      :       Indicate a specific step to run (see below)
        -h|--help                   :       Display this help message

        This script will download all 4 datasets by default. To download only
        a specific dataset, designate a step number when running the script. 
        The following step numbers are valid options:

        step1 : downloads the Human Genome Dating atlas dataset from https://human.genome.dating/
        step2 : downloads recombination maps for GRCh37 from 
                https://github.com/joepickrell/1000-genomes-genetic-maps
        step3 : downloads introgression map files from http://dical-admix.sourceforge.net/
        step4 : downloads population-specific annotation file from 
                http://ftp.ensembl.org/pub/grch37/release-103/variation/gvf/homo_sapiens/
        step5 : download GRCh37 reference genome files from 
                https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/
    "
}

if [ -z "$*" ]; then usage ; exit 1 ; fi


# Parsing command-line arguments:

OPTIONS=`getopt -o o:h --long outputdir:,runonlystep:,help -n 'download error' -- "$@"`
if [ $? != 0 ]; then echo " " ; echo "Could not parse options (see above) ..." >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTIONS"

while true ; do
    case "$1" in
        -o|--outputdir)
            case "$2" in 
                -*) echo "Please designate an output directory when using -o"; usage; exit 1 ;;
                *) outputdir=$2 ; shift 2 ;;
            esac ;;

        --runonlystep)
            case "$2" in
                -*) echo "Please designate the code4Rice3K step to execute"; usage; exit 1 ;;
                *) runonlystep=$2 ; shift 2 ;;
            esac ;;

        -h|--help)
            usage; exit 1 ;;

        --) shift ; break ;;

        *) echo "Unknown option or error" ; usage ; exit 1 ;;

    esac
done


# Control which step to run:

if [[ $runonlystep != '' ]] ; then
    steps2run=($runonlystep)
fi 

# Step 0
#
echo ""
echo "=============================================================================="
echo " Downloading datasets for the Wang etal., 2022 paper ..."
echo "=============================================================================="
echo ""

runstep1=0
runstep2=0
runstep3=0
runstep4=0
runstep5=0
for step in ${steps2run[@]} ; do
    if [[ $step == "step1" ]] ; then runstep1=1 ; fi
    if [[ $step == "step2" ]] ; then runstep2=1 ; fi
    if [[ $step == "step3" ]] ; then runstep3=1 ; fi
    if [[ $step == "step4" ]] ; then runstep4=1 ; fi
    if [[ $step == "step5" ]] ; then runstep5=1 ; fi
done


if [ ! -d "$outputdir" ] ; then
        echo ""
        echo "... new output directory $outputdir will be created in your working directory."
    	mkdir $outputdir
fi

cd $outputdir
mkdir -p Atlas Recombination Introgression Populations Reference


# Step 1
if [[ $runstep1 == 1 ]]; then
        echo ""
        echo "=============================================================================="
        echo "Step 1 Downlowding atlas SNP dating files ..."
        echo "=============================================================================="
        echo ""

	cd ${workingdir}
	cd ${outputdir}/Atlas
	
	for i in {1..22}; do
		if [[ ! -e atlas.chr${i}.csv ]]; then
			wget https://human.genome.dating/bulk/atlas.chr${i}.csv.gz && gunzip atlas.chr${i}.csv.gz
		fi
	done

	echo "Atlas files downloaded"
	echo ""
fi


# Step 2
if [[ $runstep2 == 1 ]]; then
        echo ""
        echo "=============================================================================="
        echo "Step 2 Downlowding GRCH37 recombination maps ..."
        echo "=============================================================================="
        echo ""

	cd ${workingdir}
	cd ${outputdir}/Recombination

	if [[ ! -e HapmapII_GRCh37_RecombinationHotspots.tar.gz ]]; then
		wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots/HapmapII_GRCh37_RecombinationHotspots.tar.gz && tar -xzvf HapmapII_GRCh37_RecombinationHotspots.tar.gz
	fi

	echo "Recombination maps downloaded"
	echo ""
fi


# Step 3
if [[ $runstep3 == 1 ]]; then
        echo ""
        echo "=============================================================================="
        echo "Step 3 Downlowding Neanderthal introgression files ..."
        echo "=============================================================================="
        echo ""

	cd ${workingdir}
	cd ${outputdir}/Introgression

	echo "
		Actually, I  couldn't put a script here to download this data directly.
		The Neanderthal introgression files from Steinrucken et al., 2018
		are  hosted on a google drive. You can download it yourself following
		this link here:

		https://drive.google.com/drive/folders/175ae-y9Q9Q6FQN6kQS6iGduzeAdHQZXY
		"
	echo ""
fi


# Step 4
if [[ $runstep4 == 1 ]]; then
        echo ""
        echo "=============================================================================="
        echo "Step 4 Downlowding population-specific annotations ..."
        echo "=============================================================================="
        echo ""

	cd ${workingdir}
	cd ${outputdir}/Populations

	if [[ ! -e 1000GENOMES-phase_3.gvf ]]; then
		wget http://ftp.ensembl.org/pub/grch37/release-103/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz && gunzip 1000GENOMES-phase_3.gvf.gz
	fi

	echo "Population-specific GVF file downloaded"
	echo ""
fi


# Step 5
if [[ $runstep5 == 1 ]]; then
        echo ""
        echo "=============================================================================="
        echo "Step 5 Downlowding GRCh37 reference genome ..."
        echo "=============================================================================="
        echo ""

	cd ${workingdir}
	cd ${outputdir}/Reference

	if [[ ! -e GCF_000001405.25_GRCh37.p13_genomic.fna ]]; then
		wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz && gunzip GCF_000001405.25_GRCh37.p13_genomic.fna.gz
	fi

	echo "Referece genome file downloaded"
	echo ""
fi


echo ""
echo "=============================================================================="
echo "Files downloaded."
echo "=============================================================================="
echo ""

exit

