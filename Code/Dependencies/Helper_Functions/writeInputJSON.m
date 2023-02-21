function writeInputJSON(file,current_file_contents)
    fileID = fopen(file,"w");
    fwrite(fileID,"{"+newline);
    
    current_file_fieldnames = string(fieldnames(current_file_contents));
    for name = current_file_fieldnames'
        fwrite(fileID,string(char(9))+'"'+name+'":');
        if numel(current_file_contents.(name))==1
            fwrite(fileID,num2str(current_file_contents.(name)));
        elseif numel(current_file_contents.(name))==2
            values = current_file_contents.(name);
            fwrite(fileID,"["+num2str(values(1),3)+","+num2str(values(2),3)+"]");
        elseif numel(current_file_contents.(name))==4
            values = current_file_contents.(name);
            fwrite(fileID,"[["+num2str(values(1),3)+","+num2str(values(2),3)+"],["+num2str(values(3),3)+","+num2str(values(4),3)+"]]");
        elseif numel(current_file_contents.(name))==6
            values = current_file_contents.(name);
            fwrite(fileID,"[["+num2str(values(1),3)+","+num2str(values(2),3)+","+num2str(values(3),3)+"],["+num2str(values(4),3)+","+num2str(values(5),3)+","+num2str(values(6),3)+"]]");
        end

        if ~strcmp(name,current_file_fieldnames(end))
            fwrite(fileID,","+newline);
        else
            fwrite(fileID,newline+"}");
        end
    end
    fclose(fileID);
end