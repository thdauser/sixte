()=printf("<?xml version=\"1.0\"?>\n<detector type=\"eroframestore\">\n\n<dimensions xwidth=\"384\" ywidth=\"384\"/>\n\n");
()=printf("<wcs xrpix=\"192.5\" yrpix=\"192.5\" xrval=\"0.\" yrval=\"0.\" xdelt=\"75.e-6\" ydelt=\"75.e-6\"/>\n\n");
()=printf("<response filename=\"/home/schmid/erosita/rsp/erosita_iv_7telonaxis_ff.rsp\"/>\n");
()=printf("<eventfile template=\"geneventfile.tpl\"/>\n\n");
()=printf("<threshold type=\"lower\" unit=\"keV\" value=\"150.e-3\"/>\n");
()=printf("<threshold type=\"upper\" unit=\"keV\" value=\"12.\"/>\n\n");

()=printf("<readout mode=\"time\">\n");
()=printf("\t<wait time=\"0.05\"/>\n\n");

variable i;
for(i=0; i<384; i++) {
  ()=printf("\t<readoutline lineindex=\"0\" readoutindex=\"" + string(i) + "\"/><lineshift/><wait time=\"1.e-6\"/>\n");
}

()=printf("</readout>\n");
()=printf("</detector>\n");

exit;
