function setDesktopStatus(statusText)
  % First, get the desktop reference
  try
    desktop = com.mathworks.mde.desk.MLDesktop.getInstance;     % Matlab 7+
  catch
    desktop = com.mathworks.ide.desktop.MLDesktop.getMLDesktop; % Matlab 6
  end
 
  if desktop.hasMainFrame
    % Schedule a timer to update the status text
    % Note: can't update immediately (will be overridden by Matlab's 'busy' message)
    try
      timerFcn = {@setText,desktop,statusText};
      t = timer('TimerFcn',timerFcn, 'StartDelay',0.05, 'ExecutionMode','singleShot');
      start(t);
    catch
      % Probably an old Matlab version that still doesn't have timer
      desktop.setStatusText(statusText);
    end
  else
    disp(statusText);
  end
 
%% Utility function used as setDesktopStatus's internal timer's callback
function setText(varargin)
  if nargin == 4  % just in case...
    targetObj  = varargin{3};
    statusText = varargin{4};
    targetObj.setStatusText(statusText);
  else
    % should never happen...
  end