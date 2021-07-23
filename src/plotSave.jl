#Save plot to outFile.ext for ext in exts
#(often exts=["png","pdf"])
function plotSave(p,outFile::String,exts::Array{String,1}=["png","pdf"])

  #override exts if outFile includes an extension
  for ext in [".png",".pdf",".jpg"]
    if endswith(outFile,ext)
      exts = [ chop(ext,head=1) ];
      outFile = chop(outFile,tail=length(ext)); #remove extension (will be added below)
    end
  end

  #save
  for ext in exts
    oFl = outFile*"."*ext;
    savefig(p,oFl);
    println("Wrote: $oFl");
  end
end

function plotSave(p,outFile::String,exts::String)
  plotSave(p,outFile,[exts]);
end
