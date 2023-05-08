#!/bin/sh
# Builds the image with current UID/GID, then runs the container and generates a report at mount path
# Usage: Pass (absolute) mount path and file name to the script, make sure to use bash and not sh on command line
# To be called from MTB-Report base directory - for calls from src add "cd .." before the building process

docker rm mtbreport
docker build -t mtbreport -f ./Dockerfile .

echo 'UID/GID:' $(id -u ${USER})':'$(id -g ${USER})
echo 'Using mount path:' "$1"
echo 'Input file name:' "$2"
echo '(Please make sure to pass the input file before other optional parameters)'

if [ $# -eq 0 ]; then
	echo 'Error - No parameters given. Starting the Shiny app.'
	#docker run -it -p 3838:3838 --user $(id -u ${USER}):$(id -g ${USER}) --name mtbreport mtbreport
	docker run -it -p 3838:3838 --user shiny --name mtbreport mtbreport
elif [ $# -eq 1 ]; then
	echo 'Error - Missing parameter! Make sure to pass both the mount path and at least one file as arguments.'
	#docker run -it -p 3838:3838 -v ${1}:/work/files/databases/${INPUT_FILE} --user $(id -u ${USER}):$(id -g ${USER}) --name imtb imtb /work/files/input_data/${1}

    path=${1}
	fullname="${path##*/}"
	dirname="${path%/*}"
	#basename="${fullname%.*}"
	basename="${fullname%*}"
	extension="${fullname##*.}"

	# If the file is in the same directory with the script,
	# path likely will not include any directory separator.
	#if [ "$dirname" == "$path" ]; then
  	#	dirname="."
	#fi

	# If the file has no extension, correct the variable accordingly.
	#if [ "$extension" == "$basename" ]; then
  	#	extension=""
	#fi

	#INPUT_FILE="$(basename ${1})"
	INPUT_FILE=${basename}
	echo "bind path ${1} , ${INPUT_FILE}, dirname ${dirname}"
	#docker run -it -p 3839:3838 -v ${1}:/work/files/${INPUT_FILE} --user $(id -u ${USER}):$(id -g ${USER}) --name mtbreport --rm mtbreport /work/files/${INPUT_FILE}
	docker run -it -p 3839:3838 -v ${dirname}:/work/files/ --user $(id -u ${USER}):$(id -g ${USER}) --name mtbreport --rm mtbreport /work/files/${INPUT_FILE}
elif [ $# -eq 2 ]; then
	echo "2 parameters passed"
	INPUT_DIR="$(basename ${1})"
	META_FILE="$(basename ${2})"
	docker run -it -d -p 3838:3838 -v ${1}:/work/files/ --name mtbreport mtbreport /work/files/${INPUT_FILE}
elif [ $# -eq 3 ]; then
    echo "3 parameters passed"
	docker run -it -d -p 3838:3838 -v ${1}:/work/files/ --name mtbreport mtbreport /work/files/${INPUT_FILE}
else
  echo 'Error - Too many parameters!'
fi

#echo 'Checking permissions of generated report:'
#arrIN=(${2//./ }) #Splits mount path into filename and -extension
#if [[ "$1" == */ ]] #Checks for trailing slash in mount path
#then
#	ls -l "$1"MTB_Reports/"${arrIN[0]}_report.${arrIN[1]}"
#else
#	ls -l "$1"/MTB_Reports/"${arrIN[0]}_report.${arrIN[1]}"
#fi
