variable filename = "wfi.iaat.xml";
variable file = fopen(filename, "w+");

()=fprintf(file, "<?xml version=\"1.0\"?>\n<detector type=\"wfi\" mode=\"fullframe\">\n\n");
()=fprintf(file, "<dimensions xwidth=\"64\" ywidth=\"64\"/>\n\n");
()=fprintf(file, "<wcs xrpix=\"32.5\" yrpix=\"32.5\" xrval=\"0.\" yrval=\"0.\" xdelt=\"75.e-6\" ydelt=\"75.e-6\"/>\n\n");
()=fprintf(file, "<response filename=\"/home/schmid/ixo/wfi/response/labor_new_nocompt_6enc.rsp\"/>\n");
()=fprintf(file, "<split type=\"GAUSS\" par1=\"11.e-6\"/>\n\n");
()=fprintf(file, "<eventfile template=\"geneventfile.tpl\"/>\n\n");
()=fprintf(file, "<threshold_readout_lo_keV value=\"150.e-3\"/>\n");

()=fprintf(file, "<readout mode=\"time\">\n");

variable ii;
for(ii=0; ii<64; ii++) {
  ()=fprintf(file, "\t<wait time=\"33.825e-6\"/>\n");
  ()=fprintf(file, "\t<readoutline lineindex=\"" + string(ii) + "\" readoutindex=\"" + string(ii) + "\"/>\n");
}

()=fprintf(file, "</readout>\n");
()=fprintf(file, "</detector>\n");

fclose(file);

exit;
