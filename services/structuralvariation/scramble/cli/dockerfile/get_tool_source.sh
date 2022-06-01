# Get TOOL SOURCE
echo "#[INFO] TOOL source" && \
mkdir -p $(dirname $TOOL_SOURCE) && \
if [ -e $TOOL_SOURCE ]; then \
	echo "#[INFO] TOOL TARBALL already in $TOOL_SOURCE"; \
elif $(wget --no-cache --progress=bar:force "$TOOL_SOURCE_REPO" -O $TOOL_SOURCE); then \
	echo "#[INFO] TOOL TARBALL downloaded from STARK REPO $TOOL_SOURCE_REPO"; \
	if $(wget --no-cache --progress=bar:force -nv --quiet "$(dirname $TOOL_SOURCE_REPO)/source.info" -O $(dirname $TOOL_SOURCE)/source.info); then \
		echo "#[INFO] TOOL TARBALL external source information downloaded from STARK REPO $TOOL_SOURCE_REPO " ; \
	fi ; \
elif $(wget --no-cache --progress=bar:force "$TOOL_SOURCE_EXTERNAL" -O $TOOL_SOURCE); then \
	echo "#[INFO] TOOL TARBALL downloaded from EXTERNAL SOURCE $TOOL_SOURCE_EXTERNAL"; \
	echo "$TOOL_SOURCE_EXTERNAL" > $(dirname $TOOL_SOURCE)/source.info; \
else \
	echo "#[ERROR] TOOL TARBALL NOT FOUND"; \
	exit 1; \
fi && \
if [ -e $(dirname $TOOL_SOURCE)/source.info ]; then \
	echo "#[INFO] TOOL TARBALL external source: "$(cat $(dirname $TOOL_SOURCE)/source.info) ; \
fi && \
exit 0;
