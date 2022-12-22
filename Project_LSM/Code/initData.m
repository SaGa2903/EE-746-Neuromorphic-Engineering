DATA = struct();
basepath = 'D:\PranavaStudy\IITB\sem5\EE746-Neuro\Project\Stage1\ti46\ti20\train';
count = 1;
for token = 1:10
   for speaker = 1:16
      for digit = 1:10
         relpath = '';
         if (speaker > 8)
            spknum = speaker - 8;
            spkstr = append('m',string(spknum));
            % male
         else
            spknum = speaker;
            spkstr = append('f',string(spknum));
            % female
         end
         relpath = append(relpath, '\', spkstr);
         relpath = append(relpath, '\converted');
         digitval = digit - 1;
         relpath = append(relpath, '\0', string(digitval), spkstr, 'set', string(token - 1), '.sph.wav');
         fullpath = append(basepath, relpath);
         
         if exist(fullpath)
             [y, fs] = audioread([fullpath]);
             DATA(count).sig = y/(max(abs(y)));
             DATA(count).class = digit - 1;
             count = count + 1;
         end
      end
   end
end

DATASET = 'TI46-IF-full.mat';
save(DATASET, 'DATA')