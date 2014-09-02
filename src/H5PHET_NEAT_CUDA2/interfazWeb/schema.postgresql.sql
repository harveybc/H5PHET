CREATE TABLE Admins (
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    username VARCHAR(128) NOT NULL,
    password VARCHAR(128) NOT NULL,
    email VARCHAR(128) NOT NULL
);

CREATE TABLE SNN2genom (
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    rhost INTEGER NOT NULL,
    fecha DATE NOT NULL, //FALTA EL NOW()
    indexPob INTEGER NOT NULL,
    numConex INTEGER NOT NULL,
    //TODO: FALTAN ARREGLOS DE ENABLEDS Y CONEX Y ODENEVAL
);

CREATE TABLE SNN2header (
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    rhost INTEGER NOT NULL,
    fecha DATE NOT NULL, //FALTA EL NOW()
    indexPob INTEGER NOT NULL,
    //TODO: FALTAN CAMPOS DE SNN HEADER sin arreglos verificar V2
);

CREATE TABLE Log(
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    rhost INTEGER NOT NULL,
    fecha DATE NOT NULL, //FALTA EL NOW()
    indexPobMejor INTEGER NOT NULL,
    especieMejor INTEGER NOT NULL,
    fitnessMejor INTEGER NOT NULL,
    //TODO: Faltan otros campos?
);

CREATE TABLE Config(
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    rhost INTEGER NOT NULL,
    fecha DATE NOT NULL, //FALTA EL NOW()
    sizePob INTEGER NOT NULL,
    spEspecies  INTEGER NOT NULL,

    //TODO: Faltan otros campos?
);


INSERT INTO Admins (username, password, email) VALUES ('admin', 'admin', 'harveybc@ingeni-us.com');

