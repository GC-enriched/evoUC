#!/bin/bash

while [[ $# -gt 0 ]]; do
  case $1 in
    -b|--branch)
      BRANCH="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--source)
      SOURCE="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      echo "Use -s or --source to specify the source directory with input files"
      echo "and -t or --target to specify the target directory within docker container (must be absolute path)"
      exit 1
      ;;
    --default)
      DEFAULT=YES
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
  esac
done

git push https://${ACC_TOKEN}@github.com/GC-enriched/evoUC.git