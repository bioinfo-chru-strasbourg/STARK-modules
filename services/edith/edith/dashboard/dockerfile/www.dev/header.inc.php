<?php

###############
### SESSION ###
###############

session_start();


#################
### VARIABLES ###
#################

if ($APP_SECTION=="") {
	$APP_SECTION="Dashboard";
};




#############
### LINKS ###
#############

$HEADER_LINKS='<meta charset="UTF-8">
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta name="generator" content="Mobirise v4.10.3, mobirise.com">
<meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1">
<link rel="shortcut icon" href="assets/favicon.ico" type="image/x-icon">
<meta name="description" content="'.$APP_CODE.' - '.$APP_SECTION.'">
<title>'.$APP_CODE.' - '.$APP_SECTION.'</title>
<script src="assets/plotly-latest.min.js"></script>
<link rel="stylesheet" href="assets/web/assets/mobirise-icons/mobirise-icons.css">
<link rel="stylesheet" href="assets/tether/tether.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap-grid.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap-reboot.min.css">
<link rel="stylesheet" href="assets/dropdown/css/style.css">
<link rel="stylesheet" href="assets/as-pie-progress/css/progress.min.css">
<link rel="stylesheet" href="assets/datatables/data-tables.bootstrap4.min.css">
<link rel="stylesheet" href="assets/theme/css/style.css">
<link rel="stylesheet" href="assets/gallery/style.css">
<link rel="stylesheet" href="assets/mobirise/css/mbr-additional.css" type="text/css">
<style>
	.div-wrapper {
		overflow: auto;
		max-height:400px;
		}
</style>

';



############
### MENU ###
############

$HEADER_MENU= '
<section class="menu cid-qTkzRZLJNu" once="menu" id="menu1-3">
  <nav class="navbar navbar-expand beta-menu navbar-dropdown align-items-center navbar-fixed-top collapsed bg-color transparent">
	  <button class="navbar-toggler navbar-toggler-right" type="button" data-toggle="collapse" data-target="#navbarSupportedContent" aria-controls="navbarSupportedContent" aria-expanded="false" aria-label="Toggle navigation">
		  <div class="hamburger">
			  <span></span>
			  <span></span>
			  <span></span>
			  <span></span>
		  </div>
	  </button>
	  <div class="menu-logo">
		  <div class="navbar-brand">
			  <span class="navbar-logo">
				  <a href="/">
					   <img src="assets/logo.png" alt="STARK" title="" style="height: 6rem;">
				  </a>
			  </span>
			  <span class="navbar-caption-wrap"><a class="navbar-caption text-secondary display-2" href=""> '.$APP_CODE.' - '.$APP_SECTION.'</a></span>
		 </div>
	  </div>
	  <div class="collapse navbar-collapse align-center" id="navbarSupportedContent">
	  
		  <p class="mbr-text pb-3 mbr-fonts-style display-5">
			  <img src="assets/logo.png" alt="STARK" title="" style="height: 6rem;">
		  </p>
		  <h1 class="mbr-section-title mbr-bold pb-3 mbr-fonts-style display-2">
			  '.$APP_CODE.'
		  </h1>
		  <h2 class="mbr-section-title mbr-bold pb-3 mbr-fonts-style display-5">
			  '.$APP_NAME.'
		  </h2>
		  <p class="mbr-text pb-3 mbr-fonts-style display-7">
			  Release '.$APP_RELEASE.' ['.$APP_DATE.']
			 <BR>
			  © Copyright '.$APP_COPYRIGHT.' - All Rights Reserved
			  <BR>
			  '.$APP_LICENCE.' Licence
			  <BR>
			  ['.$APP_AUTHORS.']
			  <BR>
			  <BR>
			  <a href="session.php?SESSION=unset&PHP_SELF='.$_SERVER["REQUEST_URI"].'">Reload</a>
			  <br><small><span style="color:gray">[clear data cache]</span></small>
			  <BR><BR>
		  </p>
		
	  </div>
  </nav>
</section>
';



#############
### BLANK ###
#############

$HEADER_BLANK='
<section class="header1 cid-ru7OEConn1" id="header16-1k">
	<!--
	<div class="container">
	   	<div class="row justify-content-center  align-center">
	   		<div class="card p-6 col-12 col-md-12">
	   			<div class="card-box">
	   				<h1 class="mbr-section-title mbr-bold pb-3 mbr-fonts-style display-2">
		   				<img src="assets/logo.png" width="128">
	   				</h1>
	   			</div>
	   		</div>
		</div>
	</div>
	-->
    <div class="container">
		<h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
			&nbsp;
	  	</h3>
    </div>

	</section>
';



##############
### HEADER ###
##############

$HEADER=$HEADER_LINKS.$HEADER_MENU.$HEADER_BLANK;

echo $HEADER;


?>