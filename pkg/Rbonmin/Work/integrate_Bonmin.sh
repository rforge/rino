## Usage integrate_BONMIN.sh [BONMIN.tgz]
## Integrates the latest Bonmin source package
## theussl, 2008-07-22

## where to put source files and headers
URL="http://www.coin-or.org/download/source/Bonmin/"
latest="Bonmin-0.1.4.tgz"
DESTINATION=../src/Bonmin
thirdparty_A=thirdparty_libs_A.gpg
thirdparty_B=thirdparty_libs_B.gpg

## --------------------------------------------
## Usage

usage() {
cat << "USAGE_END"
Usage: integrate_Bonmin.sh [-g ]
       integrate_Bonmin.sh [-i Bonmin_version.tgz]
       integrate_Bonmin.sh [-c]
Get, integrate or clean Bonmin sources

Options:
  -g, --get           get latest Version of Bonmin
  -i, --integrate     integrate given Bonmin sources
  -c, --clean         clean the R package's src directory

USAGE_END
        exit $1
}

## --------------------------------------------
## Read command line arguments

for x in "$@" ; do
    case "${x}" in
        -i|--integrate)        # integrate sources
             integrate=true
             ;;
        -c|--clean)            # clean sources
             clean=true
             ;;
        -g|--get)            # clean sources
             get=true
	     sources=$latest # sources gets latest
             ;;
        -t|--third_party)        # integrate sources
             third=true
             ;;
        -*)  # invalid option
             echo "$0: Invalid switch \"${x}\"!  Run $0 without parameters for help."
             exit 1
             ;;
        *)   # this should be the tarball of glpk sources
             if [[ ! -z "${sources}" ]] ; then
                 echo "$0: Only one source file allowed: \"${sources}\"!  Run $0 without parameters for help."
                 exit 1
             fi
             sources="${x}"
             ;;
    esac
done

## --------------------------------------------
## input validation

if [[ ! ( ${integrate} || ${clean} || ${get} || ${third}) ]] ; then
    echo "$0: No option given; nothing to do!"
    usage 1
    exit 1
fi

if [[ ( ${integrate} && ${clean} ) || ( ${get} && ${clean} ) || ( ${third} && ${clean} )]] ; then
    echo "$0: --clean can only be used alone!  Run $0 without parameters for help."
    exit 1
fi

if [[ -z "${sources}" && $integrate ]] ; then
    echo "$0: No source file to integrate given!"
    usage 1
    exit 1
fi

if [[ -z "${thirdparty_A}" && -z "${thirdparty_B}" && $third ]] ; then
    echo "$0: No third party libs found!"
    usage 1
    exit 1
fi


## --------------------------------------------
## integrate Bonmin sources to package

if [[ $get ]] ; then
    if [[ ! -s "${sources}" ]] ; then
	wget $URL/$sources
    else
	echo "$sources already available."
    fi
fi

if [[ $integrate ]] ; then
    
    if [[ ! -s "${sources}" ]] ; then
	echo "$0: Selected source file \"$sources\" is not available or zero!"
	usage 1
	exit 1
    fi
    Bonmin=`basename $sources .tgz`
    SOURCEDIR=${Bonmin}

    tar xzf $sources
    
    if [[ ! -d $DESTINATION ]] ; then
	mkdir $DESTINATION
    fi

    cp -r $SOURCEDIR/* $DESTINATION
    if [[ -d $SOURCEDIR ]] ; then
	rm -rf $SOURCEDIR
    fi
	
fi


if [[ $clean ]] ; then
    if [[ -d $DESTINATION ]] ; then
	rm -rf $DESTINATION
    fi
fi

if [[ $third=true ]] ; then
    if [[ -d $DESTINATION/ThirdParty/HSL ]] ; then
	gpg -o $DESTINATION/ThirdParty/HSL/ma27ad.f thirdparty_libs_A.gpg
	gpg -o $DESTINATION/ThirdParty/HSL/mc19ad.f thirdparty_libs_B.gpg

    fi
fi
             

echo "done."

exit 0
