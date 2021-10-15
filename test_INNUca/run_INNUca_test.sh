#!/usr/bin/env bash



# https://bosker.wordpress.com/2012/02/12/bash-scripters-beware-of-the-cdpath/

unset CDPATH



# https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"



CWD=$(pwd)



parse_params() {
  exports_file=''
  kraken_db=''
  threads=1
  kraken_1=''
  kraken_2=''

  while :; do
    case "${1-}" in
    --exports_file)
      exports_file="${2-}"
      shift
      ;;
    --kraken_db)
      kraken_db="${2-}"
      shift
      ;;
    --threads)
      threads="${2-}"
      shift
      ;;
    --fastq_1)
      fastq_1="${2-}"
      shift
      ;;
    --fastq_2)
      fastq_2="${2-}"
      shift
      ;;
    -?*) die "Unknown option: $1" ;;
    *) break ;;
    esac
    shift
  done

  envs_to_build=("$@")

  return 0
}

parse_params "$@"



if [[ $exports_file != '' ]]; then
  if [ -f  $exports_file ]; then
    source  $exports_file
    dependencies_option='--doNotUseProvidedSoftware'
  else
    echo "Export file was not found: $exports_file"
    exit 1
  fi
else
  # cd $DIR
  # if [ -f .exports.sh ]; then
  #   source .exports.sh
  # else
  #   echo '.exports.sh was not found'
  #   exit 1
  # fi
  # cd $CWD
  dependencies_option=''
fi



if [ -d $DIR/temp ]; then
  rm -r $DIR/temp/
fi



if [[ $fastq_1 != '' && $fastq_2 != '' ]]; then
  fastq_option="--fastq $fastq_1 $fastq_2"
else
  echo 'Downloading reads'

  mkdir -p $DIR/temp/reads

  curl -s -o $DIR/temp/reads/read_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/000/SRR5575010/SRR5575010_1.fastq.gz
  curl -s -o $DIR/temp/reads/read_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/000/SRR5575010/SRR5575010_2.fastq.gz

  echo 'Download done'

  fastq_option="--fastq $DIR/temp/reads/read_1.fastq.gz $DIR/temp/reads/read_2.fastq.gz"
fi



echo 'Running INNUca test'

if [[ $kraken_db != '' ]]; then
  kraken_options="--runKraken --krakenDB $kraken_db --krakenQuick"
else
  kraken_options=""
fi

python2 $DIR/../INNUca.py \
  --speciesExpected "Streptococcus agalactiae" \
  --genomeSizeExpectedMb 2.1 \
  $fastq_option \
  --outdir $DIR/temp/out/ \
  --threads $threads \
  $dependencies_option \
  $kraken_options \
  --runInsertSize

if [ $? -gt 0 ]; then
  echo ''
  echo 'Problem running INNUca'
  exit 1
fi

echo ''
echo 'Test INNUca done'



rm -r $DIR/temp/reads/
rm -r $DIR/temp/out/
rm -r $DIR/temp/



echo 'Everything went fine'
