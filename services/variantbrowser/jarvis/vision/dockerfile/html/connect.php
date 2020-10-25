<?php
	
	# DB1   
	/*
	$myuser='ircngs1_admin'; 
	$mypass='Adminirc1+2017!'; 
	$mydb = 'ircngs1'; 
	$myconfig= 'config.ini';
	*/	
	/*
	# DB2
	$myuser='ircngs2_admin'; 
	$mypass='Adminirc2+2017!'; 
	$mydb = 'ircngs2'; 
	$myconfig= 'config.db2.ini';
	*/
	
	# DB3
	$myuser='ircngs3_admin'; 
	$mypass='Adminirc3+2017!'; 
	# $mypass='admin3irc'; 
	$mydb = 'ircngs3'; 
	$myconfig= 'config.db3.ini';
	

	$myAnnotationConfig= 'config.annotation.ini';


	$host="hux187";

   $mydbh = mysql_connect ($host, $myuser, $mypass);
     #or die ("<h1>Could not connect to database: please try again later</h1>");
   mysql_select_db ($mydb,$mydbh);
     #or die ("<h1>Could not select proper database: please try again later</h1>");
   #echo "<h1>got connection</h1>\n"; 

# FIND information about users...



#echo USER;
#echo USERNAME;
$query="SELECT user.id AS user_id, user.username AS user_username, `group`.id AS group_id, `group`.ref AS group_ref, `project`.id AS project_id, `project`.ref AS project_ref
	FROM user
	INNER JOIN user_group ON (user.id=user_group.user_id)
	INNER JOIN `group` ON (user_group.group_id=group.id)
	INNER JOIN `project` ON (group.id=project.group_id)
	WHERE user.username='".USERNAME."'
	 "; 
#echo "<pre>"; print "$query"; echo "</pre>";
$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
while ($row = mysql_fetch_assoc($result)) {
	#echo "<pre>"; print_r($row); echo "</pre>";
	$user_id=$row['user_id'];
	$user_groups[$row['group_id']]=$row['group_ref'];
	$user_projects[$row['project_id']]=$row['group_ref']."|".$row['project_ref'];
};#while

#echo "<pre>$user_id</pre>";
#echo "<pre>"; print_r($user_groups); echo "</pre>";
#echo "<pre>"; print_r($user_projects); echo "</pre>";

#die();

 ?>
