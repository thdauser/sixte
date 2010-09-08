variable filename = "wfi.full.xml";
variable file = fopen(filename, "w+");

()=fprintf(file, "<?xml version=\"1.0\"?>\n<detector type=\"wfi\" mode=\"fullframe\">\n\n");
()=fprintf(file, "<dimensions xwidth=\"1024\" ywidth=\"1024\"/>\n\n");
()=fprintf(file, "<wcs xrpix=\"512.5\" yrpix=\"512.5\" xrval=\"0.\" yrval=\"0.\" xdelt=\"100.e-6\" ydelt=\"100.e-6\"/>\n\n");
()=fprintf(file, "<response filename=\"/home/schmid/ixo/wfi/response/ixo_wfi_default_20103103.rsp\"/>\n");
()=fprintf(file, "<split type=\"GAUSS\" par1=\"11.e-6\"/>\n\n");
()=fprintf(file, "<eventfile template=\"geneventfile.tpl\"/>\n\n");
()=fprintf(file, "<lo_keV_threshold value=\"150.e-3\"/>\n");

()=fprintf(file, "<readout mode=\"time\">\n");

variable ii;
for(ii=511; ii>=0; ii--) {
  ()=fprintf(file, "\t<wait time=\"2.e-6\"/>\n");
  ()=fprintf(file, "\t<readoutline lineindex=\"" + string(ii) + "\" readoutindex=\"" + string(ii) + "\"/>\n");
  ()=fprintf(file, "\t<readoutline lineindex=\"" + string(1023-ii) + "\" readoutindex=\"" + string(1023-ii) + "\"/>\n");
}

()=fprintf(file, "</readout>\n");
()=fprintf(file, "</detector>\n");

fclose(file);

exit;
