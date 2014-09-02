set http_proxy=http://10.11.0.1:8080
# ejecuta comando de matar automail.
d:\networking\wget.exe --proxy=on -q -t 1 -T 3 -m -nd -O d:/h5phet/tmpk0.tmp http://andi.one/phppgadmin/h5phet/kill.php
d:\networking\wget.exe --proxy=on -q -t 1 -T 3 -m -nd -O d:/h5phet/tmpk1.tmp http://credi.one/phppgadmin/h5phet/kill.php
d:\networking\wget.exe --proxy=on -q -t 1 -T 3 -m -nd -O d:/h5phet/tmpk2.tmp http://10.11.0.50/phppgadmin/h5phet/kill.php
d:\networking\wget.exe --proxy=on -q -t 1 -T 3 -m -nd -O d:/h5phet/tmpk3.tmp http://10.11.0.46/phppgadmin/h5phet/kill.php
# coloca archivos
d:\networking\winscp407 /console /script=d:\h5phet\sync.script
# ejecuta comandos de iniciar automail.
d:\networking\wget.exe --proxy=on -q -t 1 -T 3 -m -nd -O d:/h5phet/tmp0.tmp http://andi.one/phppgadmin/h5phet/syncandi.php
d:\networking\wget.exe --proxy=on -q -t 1 -T 3 -m -nd -O d:/h5phet/tmp1.tmp http://credi.one/phppgadmin/h5phet/synccredi.php
d:\networking\wget.exe --proxy=on -q -t 1 -T 3 -m -nd -O d:/h5phet/tmp0.tmp http://10.11.0.50/phppgadmin/h5phet/syncomega.php
d:\networking\wget.exe --proxy=on -q -t 1 -T 3 -m -nd -O d:/h5phet/tmp1.tmp http://10.11.0.46/phppgadmin/h5phet/syncsigma.php
