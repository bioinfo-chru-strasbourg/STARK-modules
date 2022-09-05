<?php

##################
### DATABASES ####
##################




# Databases JSON

foreach (glob ($folder_databases."/".$folder_databases_subfolder_dejavu."/STARK.database" ) as $key => $stark_database) {

	# Folders
	$stark_database_folder=dirname($stark_database);
	$stark_database_subfolder=basename($stark_database_folder);

	#echo "<br><br><br>$key $stark_database_folder $stark_database_subfolder<br>";

	# Read JSON file
	$stark_database_json=implode("",file($stark_database));
	$stark_database_obj=json_decode($stark_database_json);

	# Database infos
	$databases->{$stark_database_subfolder}->{"infos"}=$stark_database_obj;

	# Releases
	foreach (glob ($stark_database_folder."/*/STARK.database.release" ) as $key => $stark_database_release) {

		# Folders
		$stark_database_release_folder=dirname($stark_database_release);
		$stark_database_release_subfolder=basename($stark_database_release_folder);

		#echo "<br><br>$key $stark_database_release $stark_database_release_folder $stark_database_release_subfolder<br>";
		#print_r(file($stark_database_release));
		# readfile ( string $filename [, bool $use_include_path = FALSE [, resource $context ]] ) 

		# Read JSON file
		$stark_database_release_json=implode("",file($stark_database_release));
		$stark_database_release_obj=json_decode($stark_database_release_json);

		# Database infos
		#echo "<br>$stark_database_release_json "; print_r($stark_database_release_obj);
		$databases->{$stark_database_subfolder}->{"releases"}->{$stark_database_release_subfolder}=$stark_database_release_obj;

		
	};


};


// echo "<pre>";
// print_r($databases);
// echo "</pre>";


?>