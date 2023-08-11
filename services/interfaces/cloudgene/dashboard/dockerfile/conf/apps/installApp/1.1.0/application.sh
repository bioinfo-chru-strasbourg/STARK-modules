#!/bin/bash

# INPUT
ZIP=$1
echo "# Cloudgene Application ZIP file is "$(basename $ZIP)"."

# INSTALLATION
echo "# Cloudgene Application installation with ZIP file..."
if ((1)) || ! cloudgene install $ZIP; then
    echo "# Cloudgene Application installation with ZIP file FAILED."
    echo "# Cloudgene Application installation by extracting ZIP files..."
    TMP=/tmp/$RANDOM
    mkdir -p $TMP
    unzip -q $ZIP -d $TMP/
    YAML=$(ls $TMP/*/*yaml | head -n1)
    APP_ID=$(grep ^id: $YAML | head -n1 | cut -d" " -f2)
    APP_VERSION=$(grep ^version: $YAML | head -n1 | cut -d" " -f2)
    if [ $APP_ID == "" ] || [ $APP_VERSION == "" ]; then
        echo "# Cloudgene Application APP_ID or APP_VERSION DO NOT exist. Please check you YAML file."
        rm -rf $TMP
        exit 1
    else
        echo "# Cloudgene Application $APP_ID@$APP_VERSION installation..."
        APP_FOLDER="../../$APP_ID/$APP_VERSION"
        APP_FOLDER_APPS="apps/$APP_ID/$APP_VERSION"
        if [ -d $APP_FOLDER ]; then
            echo "# Cloudgene Application $APP_ID@$APP_VERSION folder exists (no change)."
        else
            echo "# Cloudgene Application $APP_ID@$APP_VERSION folder DOES NOT exist (creation)."
            mkdir -p $APP_FOLDER
            if [ -d "$TMP/$APP_VERSION" ] && (($(ls $TMP/$APP_VERSION/*yaml | wc -w))); then
                echo "# Cloudgene Application $APP_ID@$APP_VERSION folder exist in ZIP file (in $APP_VERSION subfolder)."
                cp -r $TMP/$APP_VERSION/* $APP_FOLDER/
            elif [ -d "$TMP" ] && (($(ls $TMP/*yaml | wc -w))); then
                echo "# Cloudgene Application $APP_ID@$APP_VERSION folder exist in ZIP file (in root folder)."
                cp -r $TMP/* $APP_FOLDER/$APP_VERSION/
            else
                echo "# Cloudgene Application $APP_ID@$APP_VERSION folder DOES NOT exist in ZIP file. Please check Cloudgene Application ZIP file."
                rm -rf $TMP
                exit 1
            fi;
        fi;
        YAML=$(ls $APP_FOLDER/*yaml | head -n1)
        if [ $YAML != "" ]; then
            echo "# Cloudgene Application $APP_ID@$APP_VERSION YAML is "$(basename $YAML)"."
            if ! cloudgene install $APP_FOLDER_APPS/$(basename $YAML); then
                rm -rf $TMP
                exit 1
            else
                echo "The application will be available in the interface after Cloudgene server restart."
            fi;
        else
            echo "# Cloudgene Application $APP_ID#$APP_VERSION YAML DOES NOT exist."
            rm -rf $TMP
            exit 1
        fi;
    fi;
    rm -rf $TMP
fi;

