variable filename = "epicpn.burst.xml";
variable file = fopen(filename, "w+");

()=fprintf(file, "<?xml version=\"1.0\"?>\n<detector type=\"EPIC-pn\">\n\n<dimensions xwidth=\"64\" ywidth=\"200\"/>\n\n");
()=fprintf(file, "<wcs xrpix=\"32.5\" yrpix=\"190.5\" xrval=\"0.\" yrval=\"0.\" xdelt=\"150.e-6\" ydelt=\"150.e-6\"/>\n\n");
()=fprintf(file, "<cte value=\"1\"/>\n\n");
()=fprintf(file, "<!-- TODO -->\n");
()=fprintf(file, "<response filename=\"/home/schmid/erosita/rsp/erosita_iv_1telonaxis_ff.rsp\"/>\n");
()=fprintf(file, "<!-- <split type=\"gauss\" par1=\"11.e-6\"/> -->\n\n");
()=fprintf(file, "<!-- TODO -->\n");
()=fprintf(file, "<threshold_readout_lo_keV value=\"0.\"/>\n");
()=fprintf(file, "<threshold_readout_up_keV value=\"12.\"/>\n\n");
()=fprintf(file, "<threshold_event_lo_keV value=\"150.e-3\"/>\n");
()=fprintf(file, "<threshold_split_lo_fraction value=\"0.01\"/>\n\n");
()=fprintf(file, "<eventfile template=\"geneventfile.tpl\"/>\n\n");

()=fprintf(file, "<readout mode=\"time\">\n");

variable ii;
for(ii=0; ii<200; ii++) {
  ()=fprintf(file, "\t");
  ()=fprintf(file, "<clearline lineindex=\"0\"/>");
  ()=fprintf(file, "<lineshift/>");
  ()=fprintf(file, "<wait time=\"0.72e-6\"/>");
  ()=fprintf(file, "\n");
}

for(ii=0; ii<180; ii++) {
  ()=fprintf(file, "\t<!-- TODO Calculate the event time according to Kuster (1999) -->\n");
  ()=fprintf(file, "\t");
  ()=fprintf(file, "<readoutline lineindex=\"0\" readoutindex=\"" + string(ii) + "\"/>");
  ()=fprintf(file, "<lineshift/>");
  ()=fprintf(file, "<wait time=\"0.72e-6\"/>");
  ()=fprintf(file, "\n");
}

()=fprintf(file, "</readout>\n");
()=fprintf(file, "</detector>\n");

fclose(file);

exit;
