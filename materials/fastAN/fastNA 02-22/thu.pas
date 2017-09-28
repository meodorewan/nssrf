
var arguments : integer;

begin
    if ParamCount = 0 then
    begin
         writeln( 'No parameters supplied' );
         halt(1)
    end
    else begin
         writeln('There are ', ParamCount, ' parameters' );
         for arguments := 1 to ParamCount do
             Writeln( 'Parameter ',arguments,' = ',ParamStr(arguments) );
    end
end.
