#!/bin/sh

Rscript --vanilla ./data/checkLibrary.R

cd $(dirname $0)/..
if [ $# -eq 0 ]; then
	#echo '.libPaths(c("~/lib/R/library", .libPaths())); library(shiny); library(shinyjs); runApp("src", port = 3838, host = "0.0.0.0")' | R --vanilla
	shiny-server
elif [ $# -eq 1 ]; then
	echo "Rscript --vanilla ./src/localmaster.r "$1""
	Rscript --vanilla ./src/localmaster.r "$1"
elif [ $# -eq 2 ]; then
        echo "Rscript --vanilla ./src/localmaster.r "$1" "$2""
        Rscript --vanilla ./src/localmaster.r "$1" "$2"
else
        echo "Rscript --vanilla ./src/localmaster.r "$1" "$2" "$3""
        Rscript --vanilla ./src/localmaster.r "$1" "$2" "$3"
fi
