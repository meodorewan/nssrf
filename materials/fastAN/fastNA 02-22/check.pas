{$m 100000000}
const
                max             =10000;

type            tenArr          =array[1..max,0..1000]of string[20];
                slArr           =array[1..MAX]of longint;

var             g1              :array[1..1000000,1..2]of longint;
                g2              :array[1..max,1..max]of longint;
                similar         :array[1..max,1..max]of real;
                d1,d2,r,dx      :array[1..max]of longint;
                ten1, ten2      :tenArr;
                sl1, sl2        :slArr;
                n1,n2,m1,m2     :longint;
                fi1             :string ='celeg';
                fi2             :string ='scere';


procedure docDoThi(fi:string; var ten:tenArr; kind:longint);
var f                 :text;
    s1,s2             :string;
    node,node1,node2  :integer;
 begin
   writeln('doc ',fi);

   assign(f, fi); reset(f);

   while not eof(f) do begin
       readln(f,s1);
       if pos('id',s1)>0 then begin

          if kind = 1 then inc(n1)
          else inc(n2);

          s1:=copy(s1,pos('id',s1)+3,length(s1)-pos('id',s1)-2);
          val(s1,node);
          node:=node + 1;
          readln(f,s2);
          ten[node][0]:=copy(s2,pos('label',s2)+7,length(s2)-2-pos('label',s2)-5);
          //writeln(node,' ',ten[node][0]);
       end;

       if pos('edge',s1)>0 then begin
             readln(f,s1);
             s2:=copy(s1,pos('source',s1)+6,length(s1)-pos('source',s1)-5);
             val(s2,node1);
             node1:=node1 + 1;

             readln(f,s1);
             s2:=copy(s1,pos('target',s1)+6,length(s1)-pos('target',s1)-5);
             val(s2,node2);
             node2:=node2 + 1;

             if kind = 1 then begin
                inc(m1);
                g1[m1,1]:=node1;
                g1[m1,2]:=node2;
             end
             else

             begin
                inc(m2);
                g2[node1,node2]:=1;
                g2[node2,node1]:=1;
             end;
       end;
   end;

   close(f);
 end;

procedure docSimilarSeq;
var f,g         :text;
    i,j         :longint;
    nf1,nf2     :string;
 begin
   nf1:='data\'+copy(fi1,1,2)+'-'+copy(fi2,1,2)+'\'+fi1+'-'+fi2+'.evals.pin';
   nf2:='data\'+copy(fi1,1,2)+'-'+copy(fi2,1,2)+'\'+fi2+'-'+fi1+'.evals.pin';
   writeln('doc ',nf1);

   assign(f,nf1); reset(f);
   assign(g,nf2); rewrite(g);
   while not seekeof(f) do begin
      readln(f,i,j,similar[i+1][j+1]);
      writeln(g,j,' ',i,' ',similar[i+1][j+1]);
   end;
   close(f);
   close(g);
 end;

function layID(sname :string; var ten:tenArr; n:longint):longint;
var i :longint;
 begin
   for i:=1 to n do
     if ten[i][0] = sname then exit(i);
   exit(0);
 end;

procedure layGO(fi:string; var ten:tenArr; var sl:slArr; n:longint);
var g       :text;
    s1,s2   :ansistring;
    node,v  :integer;
 begin
   writeln('doc ',fi);
   assign(g, fi);
   reset(g);

   while not eof(g) do begin

      readln(g,s1);

      if pos('GO',s1) = 0 then continue;

      s2:=copy(s1,1,pos('|',s1)-1);
      node:=layID(s2, ten, n);

      if node=0 then continue;

      //write(node,'(',s2,')','|');

      while (pos('GO',s1)) > 0 do begin
         v:=pos('GO',s1);
         delete(s1,1,v-1);

         if pos('|',s1)>0 then s2:=copy(s1,1,pos('|',s1)-1)
         else s2:=s1;

         inc(sl[node]);
         ten[node][sl[node]]:=s2;
         delete(s1,1,length(s2));
         //write(s2,'|');
      end;
      //writeln;

   end;

   close(g);

 end;

function hop(node1, node2:longint):longint;
var i,j,ok,dem:longint;
 begin
   dem:=sl2[node2];
   for i:=1 to sl1[node1] do
     begin
       ok:=1;
       for j:=1 to sl2[node2] do
         if ten1[node1][i] = ten2[node2][j] then ok:=0;
       dem:=dem + ok;
     end;
   exit(dem);
 end;

function giao(node1, node2:longint):longint;
var i,j,ok,dem:longint;
 begin
   dem:=0;
   for i:=1 to sl1[node1] do
     begin
       ok:=0;
       for j:=1 to sl2[node2] do
         if ten1[node1][i] = ten2[node2][j] then ok:=1;
       dem:=dem + ok;
     end;
   exit(dem);
 end;

function tinhGOS:extended;
var i,j,k,t1,t2 :longint;
    gos         :extended;
 begin
   gos:=0;
   for k:=1 to n1 do begin
      i:=k;
      j:=r[k];
      if j=0 then continue;

      t1:=giao(i,j);
      t2:=hop(i,j);

      if t2>0 then gos:=gos + t1/t2;
   end;
   exit(gos);
 end;

function tinhGOSSPINAL(c:char):extended;
var g     :text;
    s1,s2 :string;
    gos   :extended;
    t1,t2,node1,node2:longint;
 begin
    assign(g,'data\'+copy(fi1,1,2)+'-'+copy(fi2,1,2)+'\'
              +fi1+'-'+fi2+'-spinalmatch-original'+c);
    reset(g);
    readln(g,s1);
    readln(g,s1);

    gos:=0;
    while not eof(g) do begin
       readln(g,s1);
       s2:=copy(s1,1,pos(' ',s1)-1);
       node1:=layID(s2,ten1,n1);

       //write(s2,'|',node1,'|');

       s2:=copy(s1,pos(' ',s1)+1,length(s1)-pos(' ',s1));
       node2:=layID(s2,ten2,n2);
       //writeln(s2,'|',node2,'|');

       t1:=giao(node1,node2);
       t2:=hop(node1,node2);

       if t2>0 then gos:=gos + t1/t2;

    end;

    close(g);

    exit(gos);
 end;


procedure thongKeACOGA(fi:string);
var f                   :text;
    i,j,k               :longint;
    err                        :string;
    s,alpha                    :extended;
    svung,sdinh                :longint;
    maxsdinh,maxscanh,lvung    :longint;
    m,mt                       :longint;

   procedure go(i:integer);
   var j   :longint;
   begin
     dx[i]:=svung; inc(sdinh);
     for j:=1 to n2 do
       if (dx[j]=0)and(g2[i,j]=2) then begin
          go(j);
       end;
   end;

 begin
   err:='Dung';
   s:=0;
   assign(f,fi); reset(f);
   readln(f,m);
   while not seekeof(f) do begin
      readln(f,i,j);

      if (i<0)or(i>=n1)or(d1[i+1]=1) then
        begin err:='chi so dinh sai'; break; end;

      d1[i+1]:=1;
      if j=-1 then begin r[i+1]:=-1; continue; end;

      if (j<0)or(j>=n2)or(d2[j+1]=1) then
        begin err:='chi so dinh sai'; break; end;

      r[i+1]:=j+1;

      d2[j+1]:=1;
      s:= s + similar[i+1,j+1];
   end;
   close(f);

   mt:=0;
   for k:=1 to m1 do begin
        if (r[g1[k,1]]>0) and (r[g1[k,2]]>0) then
           mt:=mt + g2[r[g1[k,1]],r[g1[k,2]]];
   end;

   writeln(err);
   writeln(m);
   writeln(mt);

   for k:=3 to 7 do begin
        alpha:=k / 10;
        writeln(alpha:0:2,': ',alpha*mt+(1-alpha)*s:0:5);
   end;

   for k:=1 to m1 do begin
        if (r[g1[k,1]]>0) and (r[g1[k,2]]>0) and (g2[r[g1[k,1]],r[g1[k,2]]]=1) then begin
                g2[r[g1[k,1]],r[g1[k,2]]]:=2;
                g2[r[g1[k,2]],r[g1[k,1]]]:=2;
        end;
   end;

   svung:=0;
   maxsdinh:=0;
   maxscanh:=0;

   for k:=1 to n2 do
     if dx[k]=0 then begin
        inc(svung);
        sdinh:=0;

        go(k);

        if sdinh>maxsdinh then begin
                maxsdinh:=sdinh;
                lvung:=svung;
        end;
     end;

   writeln('so vung : ',svung);
   writeln('so nut : ',maxsdinh);
   for i:=1 to n2 do
     if dx[i]=lvung then
     for j:=i+1 to n2 do
       if (dx[j]=lvung)and(g2[i,j]=2) then inc(maxscanh);
   writeln('so canh : ',maxscanh);

   writeln;

   writeln('ACOGA -> ',tinhGOS:0:5);
   writeln;
 end;

BEGIN
   if ParamCount > 0 then
    begin
        fi1:=paramStr(1);
        fi2:=paramStr(2);

        if pos('.ga',fi1)>0 then fi1:=copy(fi1,1,pos('.ga',fi1)-1);
        if pos('.ga',fi2)>0 then fi2:=copy(fi2,1,pos('.ga',fi2)-1);

    end;

   docDoThi('data\'+fi1+'\'+fi1+'.tab.gml', ten1, 1);
   writeln(n1,' ',m1);
   docDoThi('data\'+fi2+'\'+fi2+'.tab.gml', ten2, 2);
   writeln(n2,' ',m2);

   docSimilarSeq;

   layGO('data\'+fi1+'\'+fi1+'_gene_to_go.txt', ten1, sl1, n1);
   layGO('data\'+fi2+'\'+fi2+'_gene_to_go.txt', ten2, sl2, n2);

   thongkeACOGA(fi1+'.ga_'+fi2+'.ga.out');

   writeln('SPINAL 0.3 -> ',tinhGOSSPINAL('3'):0:5);
   writeln('SPINAL 0.4 -> ',tinhGOSSPINAL('4'):0:5);
   writeln('SPINAL 0.5 -> ',tinhGOSSPINAL('5'):0:5);
   writeln('SPINAL 0.6 -> ',tinhGOSSPINAL('6'):0:5);
   writeln('SPINAL 0.7 -> ',tinhGOSSPINAL('7'):0:5);
END.
