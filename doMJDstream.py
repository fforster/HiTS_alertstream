import os, re, sys
import matplotlib.pyplot as plt
import numpy as np
import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter
from astropy.io import fits

alerts = {}

reader = DataFileReader(open("hits-demo.avro", "rb"), DatumReader())
for alert in reader:
    MJD = alert["sciMJD"]
    print(MJD, alert["Id"])
    if not MJD in alerts.keys():
        alerts[MJD] = []
    alerts[MJD].append(alert)
reader.close()

schema = avro.schema.parse(open("hits.avsc", "rb").read())

for MJD in sorted(alerts.keys()):
    writer = DataFileWriter(open("alerts/hits-demo_%s.avro" % MJD, "wb"), DatumWriter(), schema)
    for alert in alerts[MJD]:
        writer.append(alert)
    writer.close()
