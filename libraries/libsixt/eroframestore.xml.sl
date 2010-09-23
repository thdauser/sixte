variable filename = "eroframestore.xml";
variable file = fopen(filename, "w+");

()=fprintf(file, "<?xml version=\"1.0\"?>\n");
()=fprintf(file, "<instrument type=\"eROSITA\">\n\n");

()=fprintf(file, "<telescope>\n");
()=fprintf(file, "<focallength value=\"1.6\"/>\n");
()=fprintf(file, "<fov diameter=\"1.02\"/>\n");
()=fprintf(file, "<psf filename=\"/home/schmid/erosita/psf/erosita.psf.7.5mum.fits\"/>\n");
()=fprintf(file, "</telescope>\n\n");

()=fprintf(file, "<detector type=\"framestoreCCD\">\n");
()=fprintf(file, "<dimensions xwidth=\"384\" ywidth=\"384\"/>\n");
()=fprintf(file, "<wcs xrpix=\"192.5\" yrpix=\"192.5\" xrval=\"0.\" yrval=\"0.\" xdelt=\"75.e-6\" ydelt=\"75.e-6\"/>\n");
()=fprintf(file, "<cte value=\"1\"/>\n");
()=fprintf(file, "<response filename=\"/home/schmid/erosita/rsp/erosita_iv_1telonaxis_ff.rsp\"/>\n");
()=fprintf(file, "<split type=\"exponential\" par1=\"0.355\"/>\n");
()=fprintf(file, "<threshold_readout_lo_keV value=\"0.\"/>\n");
()=fprintf(file, "<threshold_readout_up_keV value=\"12.\"/>\n");
()=fprintf(file, "<threshold_event_lo_keV value=\"150.e-3\"/>\n");
()=fprintf(file, "<threshold_split_lo_fraction value=\"0.01\"/>\n");
()=fprintf(file, "<eventfile template=\"geneventfile.tpl\"/>\n\n");

()=fprintf(file, "<readout mode=\"time\">\n");
()=fprintf(file, "\t<wait time=\"0.05\"/>\n\n");

variable ii;
for(ii=0; ii<384; ii++) {
  ()=fprintf(file, "\t<readoutline lineindex=\"0\" readoutindex=\"" + string(ii) + "\"/><lineshift/>");
  ()=fprintf(file, "<wait time=\"1.e-6\"/>\n");
}

()=fprintf(file, "</readout>\n");
()=fprintf(file, "</detector>\n");
()=fprintf(file, "</instrument>\n");

fclose(file);

exit;
