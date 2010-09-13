variable filename = "eroframestore.xml";
variable file = fopen(filename, "w+");

()=fprintf(file, "<?xml version=\"1.0\"?>\n<detector type=\"eroframestore\">\n\n<dimensions xwidth=\"384\" ywidth=\"384\"/>\n\n");
()=fprintf(file, "<wcs xrpix=\"192.5\" yrpix=\"192.5\" xrval=\"0.\" yrval=\"0.\" xdelt=\"75.e-6\" ydelt=\"75.e-6\"/>\n\n");
()=fprintf(file, "<cte value=\"1\"/>\n\n");
()=fprintf(file, "<response filename=\"/home/schmid/erosita/rsp/erosita_iv_1telonaxis_ff.rsp\"/>\n");
()=fprintf(file, "<split type=\"exponential\" par1=\"0.355\"/>\n\n");
()=fprintf(file, "<lo_keV_threshold value=\"150.e-3\"/>\n");
()=fprintf(file, "<up_keV_threshold value=\"12.\"/>\n\n");
()=fprintf(file, "<eventfile template=\"geneventfile.tpl\"/>\n\n");

()=fprintf(file, "<readout mode=\"time\">\n");
()=fprintf(file, "\t<wait time=\"0.05\"/>\n\n");

variable ii;
for(ii=0; ii<384; ii++) {
  ()=fprintf(file, "\t<readoutline lineindex=\"0\" readoutindex=\"" + string(ii) + "\"/><lineshift/>\n");
  ()=fprintf(file, "<wait time=\"1.e-6\"/>\n");
}

()=fprintf(file, "</readout>\n");
()=fprintf(file, "</detector>\n");

fclose(file);

exit;
