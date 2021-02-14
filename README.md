##  Respuesta a la infección por SARS-CoV2
En este trabajo estudiaremos la respuesta desarrollada en las células del epitelio del pulmón a la infección por SARS-Cov2 mediante el análisis de los perfiles de expresión génica publicados en el dataset GEO GSE147507. Estos perfiles de expresión se analizarán mediante el modelado de redes de coexpresión génica, con el paquete de R WGCNA.

Somos estudiantes de grado de Ingeniería de la Salud con mención en Bioinformática, de la universidad de Málaga. Este proyecto pertenece a un trabajo de la asignatura de Biología de Sistemas.

Componentes:
* Ana Galiá Caravaca
* Laura Fernández García
* Mariana González Jiménez
* Juan Sánchez Rodríguez 
* Susana Cano Marín (coordinadora)

**Instalación**
Para la ejecución de este proyecto sera necesario:
- Tener instalado R con la versión 3.0.4
- Paquete de R **WGCNA**, paquete principal de R empleado en el desarrollo de este trabajo, para el análisis de redes de coexpresión de genes ponderados.

Se adjunta un fichero "WGCNA_deps_R.R" para la actualización de la versión de R necesaria, para Windows. Así como la instalación de los paquetes de R empleados.

**Documentación**
Los datos GEO GSE147507 han sido descargados NCBI y se facilitan en la carpeta de "code" (GSE147507_RawReadCounts_Human.tsv).

**Flujo de trabajo**
- setup.sh: actualiza R y carga librerías.
- launch.sh: lee setup.sh y ejecuta los ficheros de R.
	
¡¡¡En el flujo de trabajo solo están cargados los ficheros "WGCNA_deps_R.R", "1_PrecProcClust.R" y "2_AnalisisCoexWGCNA.R". Debido al tiempo de ejecución. !!!


Ejemplo de ejecución en GIT BASH:

sh launch.sh "\carpeta_de_trabajo" "\carpeta_code" "\carpeta_de_trabajo\software" "\carpeta_results"

Donde el contenido de los comillas son las rutas de las carpetas y la carpeta de las dependencias es la llamada "software" que crearemos en el directorio de trabajo automáticamente en setup.sh


