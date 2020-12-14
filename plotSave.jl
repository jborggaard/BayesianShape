#Save plot to outFile.ext for ext in exts
#(often exts=["png","pdf"])
function plotSave(p,outFile::String,exts::Array{String,1})
  for ext in exts
    oFl = outFile*"."*ext;
    savefig(p,oFl);
    println("Wrote: $oFl");
  end
end
function plotSave(p,outFile::String,exts::String)
  plotSave(p,outFile,[exts]);
end
