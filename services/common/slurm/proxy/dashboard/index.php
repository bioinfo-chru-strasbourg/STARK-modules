<?php
error_reporting(E_ALL);
session_start();

if (isset($_GET) && isset($_GET['user'])) {
	session_start();
	$_SESSION['user_name'] = $_GET['user'];
	header("X-SLURM-USER-NAME: ".$_SESSION['user_name']);
	#echo "<p>Hello {$_GET['user']}.</p>";

	if (!isset($_GET['token']) || $_GET['token'] == "") {
		unset($_SESSION['user_token']);
		#echo "<p>Using slurm user for authentication proxy.</p>";
	} else {
		$_SESSION['user_token'] = $_GET['token'];
		header("X-SLURM-USER-TOKEN: ".$_SESSION['user_token']);
		#echo "<p>You entered {$_GET['token']} as your token.</p>";
	}

} else {

	header('HTTP/1.0 401 Unauthorized');
	unset($_SESSION['user_name']);
	unset($_SESSION['user_token']);

};

#include "session.php";

?>

<html>

<head>
	<title>Dashboard</title>
	<link rel="stylesheet" type="text/css" href="main.css">
<!--
	<link rel="stylesheet" type="text/css" href="main.css">
	<link rel="stylesheet" type="text/css" href="handsontable/dist/handsontable.full.min.css">
	<link rel="stylesheet" type="text/css" href="handsontable/static/css/main.css">
	<script src="handsontable/dist/handsontable.full.min.js"></script>
-->
</head>

<body>

<?php

if (session_status() == PHP_SESSION_NONE || !isset($_SESSION['user_name']) || $_SESSION['user_name'] == "") {

?>

<!--
<form action="?" method="get">
<div class="imgcontainer">
    
  </div>

  <div class="container">
    <label for="uname"><b>Username</b></label>
    <input type="text" placeholder="Enter Username" name="uname" required>

    <label for="psw"><b>Password</b></label>
    <input type="password" placeholder="Enter Password" name="psw" required>

    <button type="submit">Login</button>
    <label>
      <input type="checkbox" checked="checked" name="remember"> Remember me
    </label>
  </div>

  <div class="container" style="background-color:#f1f1f1">
    <button type="button" class="cancelbtn">Cancel</button>
    <span class="psw">Forgot <a href="#">password?</a></span>
  </div>
</form>

-->


<form action="?" method="get">
Login/Password	<input type="text" name="user"> <input type="text" name="token"> <input type="submit" value='login'>

<span title='Authentication Options:
- Per user token:
	User: "fred"
	Password: use generated token from scontrol
- Authentication Proxy:
	User: "fred"
	Password: leave empty to use "slurm"
		user as an authentication proxy.
	'>help
</span>

</form>



<?php

die();

} else {

	$login_message = "<form action='?' method='get'> ";
	$login_message .= "Hello {$_SESSION['user_name']} ";
	if ($_SESSION['user_token']) {
		$login_message .= "(You entered {$_SESSION['user_token']} as your token) ";
	} else {
		$login_message .= "(Using slurm user for authentication proxy) ";
	};
	$login_message .= "<input type='submit' value='logout'></form>";
};

?>

<p>
<?php print($login_message); ?>
</p>


<p>
	<A href='sacct.php' target='sacct'>List of jobs</a>
</p>


<?php
//include "sacct.php";
?>


</body>
</html>
