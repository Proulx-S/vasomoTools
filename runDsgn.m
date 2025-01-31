classdef runDsgn
    % Encapsulate definition of stimulus design.

   properties
      % Header info
      task      % [str] label
      onsetList % [num] onset time of each event/trial, in seconds (trigger time = 0)
      ondurList % [num] duration of each event/trial, in seconds
      dt        % [num] time grid on which events/trials were defined (not sure it's usefull)
      cond      % [int] index of event/trial condition (0 is special for null event/trial, other conditions should be increments of 1)
      condLabel % [cellstr] event/trial condition labels (indices corresponds to values of cond)
   end
end