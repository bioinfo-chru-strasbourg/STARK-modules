
function merge_filters($filter1,$filter2,$filter_array) {
	# if filter 1 is based on filter 2, ad filter 2 not based on filter1
	#if (isset($filter_array[$filter1]["basedon"]) && $filter_array[$filter2]["basedon"]!=$filter1) {
		# if filter 2 is based on another filter
		#echo "$filter1,$filter2,".$filter_array[$filter1]["basedon"]."<BR>";
		if (isset($filter_array[$filter2]["basedon"]) && isset($filter_array[$filter_array[$filter2]["basedon"]]) && $filter_array[$filter2]["basedon"]!=$filter1) {
			$filter_array=merge_filters($filter2,$filter_array[$filter2]["basedon"],$filter_array);
		};#if
		if (isset($filter_array[$filter1]["basedon"]) && isset($filter_array[$filter_array[$filter1]["basedon"]]) && $filter_array[$filter1]["basedon"]!=$filter2) {
			$filter_array=merge_filters($filter1,$filter_array[$filter1]["basedon"],$filter_array);
		};#if
		# merge arrays
		$filter_array[$filter1]=array_merge($filter_array[$filter2],$filter_array[$filter1]);
	#};#if
	return $filter_array;
}

function basedon_filters($filter_array) {

	foreach ($filter_array as $filter => $filter_criteria) {
		#echo "FILTER: $filter<BR>";
		if (isset($filter_array[$filter]["basedon"]) && isset($filter_array[$filter_array[$filter]["basedon"]])) {
			$filter_array=merge_filters($filter,$filter_array[$filter]["basedon"],$filter_array);
		};#if
	};#foreach
	return $filter_array;

};#function merge_filters

function load_filter($filters,$filter_array,$filter_option_array=array()) {

	foreach ($filters as $filter) {
		foreach ($filter_array[$filter] as $annotation => $annotation_filter) {
			#echo "$annotation => $annotation_filter<BR>";	
			if (is_array($annotation_filter)) {
				foreach ($annotation_filter as $annotation_filter_key => $annotation_filter_value) {
					#echo "&nbsp;&nbsp;&nbsp;&nbsp;$annotation_filter_key => $annotation_filter_value<BR>";
					$filter_option_array[]="$annotation:$annotation_filter_value";
				};#foreach
			} else {
				$filter_option_array[]="$annotation:$annotation_filter";
			};#if
		};#foreach
	};#foreach
	#print_r($filter_option_array);
	return $filter_option_array;

};#function load_filter
