echo "#[INFO] TOOL cleaning" && \
rm -rf $TOOL_SOURCE_BUILD && \
if (($REMOVE_SOURCES)); then rm -rf $SOURCES/$SOURCES_FOLDER/tools/$TOOL_NAME; fi && \
echo "#[INFO] TOOL $TOOL_NAME/$TOOL_VERSION installed"