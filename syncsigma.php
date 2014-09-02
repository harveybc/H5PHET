<HTML>
<HEAD>
<TITLE>sync</TITLE>
</HEAD>
<BODY>
	Syncing.<br>
	<?php 
		passthru("/home/harveybc/h5/autoscript");
		passthru("/usr/sbin/automail -t96 -e8 -n8 -d5 -p/var/www/phppgadmin/h5phet/ -x/usr/sbin/automailconf &");
	?><br>
	Synced.<br>
</BODY>
</HTML>
