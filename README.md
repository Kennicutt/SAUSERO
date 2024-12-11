SAUSERO es una rutina de reducción para el modo imagen de OSIRIS+. Hay otra versión válida para
la reducción de las observaciones adquiridas con el anterior sensor.

Para su uso, se debe llevar a cabo las siguientes acciones:

1. Conectarse a Phase3. Esto puede hacerse en remoto a través del siguiente comando:
"ssh -XY phase3@calp-drm".

2. Activar el entorno de SAUSERO mediante el siguiente comando: "conda activate sausero".

3. Para el ejecución de la rutina se comanda: "SAUSERO <PRG> <OB>".

4. Una vez finalizada la reducción, se indicará la ruta y la lista de imágenes finales generadas.

Los resultados generados son una imagen por cada banda observada. En dicha imagen se habrá hecho
sustracción de bias, correción de flatfield (y corrección de fringing en caso Sloan z), alineado
de las imágenes entre sí y sumado, calibración astrométrica y calibración fotométrica (estimación del
ZP +- error). Para la eliminación de efectos cosméticos se hace uso de una BPM y se aplica el
algoritmo de LACosmic para los rayos cósmicos.