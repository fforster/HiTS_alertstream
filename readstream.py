import os, re, sys
import matplotlib.pyplot as plt
import numpy as np
import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter
from astropy.io import fits


showschema = True

for alertfile in sorted(os.listdir("alerts")):
    
    print(alertfile)
    reader = DataFileReader(open("alerts/%s" % alertfile, "rb"), DatumReader())
    
    # show schema once
    if showschema:
        schema = reader.datum_reader.writers_schema
        for field in schema.fields:
            print(field.name, field.doc)
        showschema = False

    # show alert information
    for alert in reader:
        print(alert["Id"])
        #print(alert.keys())


reader.close()
