/********************** 
 * Roar-Sharable Test *
 **********************/

import { PsychoJS } from './lib/core-2020.1.js';
import * as core from './lib/core-2020.1.js';
import { TrialHandler } from './lib/data-2020.1.js';
import { Scheduler } from './lib/util-2020.1.js';
import * as util from './lib/util-2020.1.js';
import * as visual from './lib/visual-2020.1.js';
import * as sound from './lib/sound-2020.1.js';

// init psychoJS:
const psychoJS = new PsychoJS({
  debug: true
});

// open window:
psychoJS.openWindow({
  fullscr: true,
  color: new util.Color([0, 0, 0]),
  units: 'height',
  waitBlanking: true
});

// store info about the experiment session:
let expName = 'ROAR-Sharable';  // from the Builder filename that created this script
let expInfo = {'participant': '', 'password': ''};

// schedule the experiment:
psychoJS.schedule(psychoJS.gui.DlgFromDict({
  dictionary: expInfo,
  title: expName
}));

const flowScheduler = new Scheduler(psychoJS);
const dialogCancelScheduler = new Scheduler(psychoJS);
psychoJS.scheduleCondition(function() { return (psychoJS.gui.dialogComponent.button === 'OK'); }, flowScheduler, dialogCancelScheduler);

// flowScheduler gets run if the participants presses OK
flowScheduler.add(updateInfo); // add timeStamp
flowScheduler.add(experimentInit);
flowScheduler.add(exSetUpRoutineBegin());
flowScheduler.add(exSetUpRoutineEachFrame());
flowScheduler.add(exSetUpRoutineEnd());
const InstLoop1LoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(InstLoop1LoopBegin, InstLoop1LoopScheduler);
flowScheduler.add(InstLoop1LoopScheduler);
flowScheduler.add(InstLoop1LoopEnd);
const pracBlockLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(pracBlockLoopBegin, pracBlockLoopScheduler);
flowScheduler.add(pracBlockLoopScheduler);
flowScheduler.add(pracBlockLoopEnd);
flowScheduler.add(clericUnlockRoutineBegin());
flowScheduler.add(clericUnlockRoutineEachFrame());
flowScheduler.add(clericUnlockRoutineEnd());
const InstLoop2LoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(InstLoop2LoopBegin, InstLoop2LoopScheduler);
flowScheduler.add(InstLoop2LoopScheduler);
flowScheduler.add(InstLoop2LoopEnd);
const threeBlocksLoopLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(threeBlocksLoopLoopBegin, threeBlocksLoopLoopScheduler);
flowScheduler.add(threeBlocksLoopLoopScheduler);
flowScheduler.add(threeBlocksLoopLoopEnd);
flowScheduler.add(gameEndRoutineBegin());
flowScheduler.add(gameEndRoutineEachFrame());
flowScheduler.add(gameEndRoutineEnd());
flowScheduler.add(CreditsRoutineBegin());
flowScheduler.add(CreditsRoutineEachFrame());
flowScheduler.add(CreditsRoutineEnd());
flowScheduler.add(quitPsychoJS, '', true);

// quit if user presses Cancel in dialog box:
dialogCancelScheduler.add(quitPsychoJS, '', false);

psychoJS.start({
  expName: expName,
  expInfo: expInfo,
  });


var frameDur;
function updateInfo() {
  expInfo['date'] = util.MonotonicClock.getDateStr();  // add a simple timestamp
  expInfo['expName'] = expName;
  expInfo['psychopyVersion'] = '2020.1.2';
  expInfo['OS'] = window.navigator.platform;

  // store frame rate of monitor if we can measure it successfully
  expInfo['frameRate'] = psychoJS.window.getActualFrameRate();
  if (typeof expInfo['frameRate'] !== 'undefined')
    frameDur = 1.0 / Math.round(expInfo['frameRate']);
  else
    frameDur = 1.0 / 60.0; // couldn't get a reliable measure so guess

  // add info from the URL:
  util.addInfoFromUrl(expInfo);
  
  return Scheduler.Event.NEXT;
}


var exSetUpClock;
var currentBlock;
var nNotLast;
var coins;
var corrAnimMov;
var ncorrectSound;
var nincorrectSound;
var nhalfway;
var InstLoopClock;
var instloop_resp;
var instImg;
var pageloop;
var LexBlockClock;
var lexCrossStim;
var word_stim;
var lexCrossStimAfter;
var lex_resp;
var lexArrowKey;
var pCorrectSoundCodeClock;
var correctSoundPracClock;
var coin_sound_2;
var triArrowKeyCorr_2;
var lexCrossStimCSound;
var incorrectPracClock;
var incorrectPrac_resp;
var incorrect_image;
var fail_soundPrac;
var correction_audio;
var clericUnlockClock;
var clericUnlock_resp;
var clericUnlockImg;
var clericUnlockMovClock;
var clericUnlockMov;
var clericUnlockWav;
var blockSetUpClock;
var score;
var blockCount;
var animMov;
var unlockMov;
var charUnlockImg;
var charUnlockWav;
var valleyWav;
var valleyImg;
var nPlayEnd;
var correctSoundCodeClock;
var correctSoundClock;
var coin_sound_2lex;
var lexARrowKeyCorr_2;
var lexCrossStimCSoundLex;
var inCorrSoundClock;
var coin_sound_3;
var arrowKeyLexCorr;
var lexCrossStim1;
var blockAnimLoopClock;
var coin_icon;
var correctAnimText;
var animCoinSound;
var fixateOnCrossClock;
var lexCross_b4;
var lexArrowKey_again;
var halfwayBreakClock;
var halfway_resp;
var halfwayImg;
var halfwayWav;
var halfwaySoundWav;
var InnInstClock;
var valleyinst_img;
var valley_resp;
var inn_voice;
var charUnlockClock;
var CharUn_resp;
var CharUnlockImg;
var amountofGold;
var charUn_voice;
var valleyInstClock;
var valleyinst_image;
var valley_resp1;
var valley_voice;
var gameEndClock;
var theEnd_text;
var endMovClock;
var endMov;
var EndImg;
var end_voice;
var key_resp_end;
var CreditsClock;
var Credits_text;
var credits_wav;
var credits_resp;
var globalClock;
var routineTimer;
function experimentInit() {
  // Initialize components for Routine "exSetUp"
  exSetUpClock = new util.Clock();
  currentBlock = 0;
  nNotLast = 0;
  coins = 0;
  corrAnimMov = "Block1.mp4";
  ncorrectSound = 0;
  nincorrectSound = 0;
  nhalfway = 0;
  
  // Initialize components for Routine "InstLoop"
  InstLoopClock = new util.Clock();
  instloop_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  instImg = new visual.ImageStim({
    win : psychoJS.window,
    name : 'instImg', units : 'pix', 
    image : undefined, mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  pageloop = new sound.Sound({
    win: psychoJS.window,
    value: 'A',
    secs: (- 1),
    });
  pageloop.setVolume(1);
  // Initialize components for Routine "LexBlock"
  LexBlockClock = new util.Clock();
  lexCrossStim = new visual.ShapeStim ({
    win: psychoJS.window, name: 'lexCrossStim', 
    vertices: 'cross', size:[0.015, 0.015],
    ori: 0, pos: [0, 0],
    lineWidth: 1, lineColor: new util.Color([1, 1, 1]),
    fillColor: new util.Color([1, 1, 1]),
    opacity: 1, depth: 0, interpolate: true,
  });
  
  word_stim = new visual.TextStim({
    win: psychoJS.window,
    name: 'word_stim',
    text: 'default text',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0], height: 0.045,  wrapWidth: undefined, ori: 0,
    color: new util.Color('white'),  opacity: 1,
    depth: -1.0 
  });
  
  lexCrossStimAfter = new visual.ShapeStim ({
    win: psychoJS.window, name: 'lexCrossStimAfter', 
    vertices: 'cross', size:[0.015, 0.015],
    ori: 0, pos: [0, 0],
    lineWidth: 1, lineColor: new util.Color([1, 1, 1]),
    fillColor: new util.Color([1, 1, 1]),
    opacity: 1, depth: -2, interpolate: true,
  });
  
  lex_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  lexArrowKey = new visual.ImageStim({
    win : psychoJS.window,
    name : 'lexArrowKey', units : undefined, 
    image : 'arrowkey_lex.png', mask : undefined,
    ori : 0, pos : [0, (- 0.2)], size : [0.2, 0.2],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -4.0 
  });
  // Initialize components for Routine "pCorrectSoundCode"
  pCorrectSoundCodeClock = new util.Clock();
  // Initialize components for Routine "correctSoundPrac"
  correctSoundPracClock = new util.Clock();
  coin_sound_2 = new sound.Sound({
    win: psychoJS.window,
    value: 'coin_sound.wav',
    secs: (- 1),
    });
  coin_sound_2.setVolume(0.5);
  triArrowKeyCorr_2 = new visual.ImageStim({
    win : psychoJS.window,
    name : 'triArrowKeyCorr_2', units : undefined, 
    image : 'arrowkey_lex.png', mask : undefined,
    ori : 0, pos : [0, (- 0.2)], size : [0.2, 0.2],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  lexCrossStimCSound = new visual.ShapeStim ({
    win: psychoJS.window, name: 'lexCrossStimCSound', 
    vertices: 'cross', size:[0.015, 0.015],
    ori: 0, pos: [0, 0],
    lineWidth: 1, lineColor: new util.Color([1, 1, 1]),
    fillColor: new util.Color([1, 1, 1]),
    opacity: 1, depth: -2, interpolate: true,
  });
  
  // Initialize components for Routine "incorrectPrac"
  incorrectPracClock = new util.Clock();
  incorrectPrac_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  incorrect_image = new visual.ImageStim({
    win : psychoJS.window,
    name : 'incorrect_image', units : 'pix', 
    image : undefined, mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  fail_soundPrac = new sound.Sound({
    win: psychoJS.window,
    value: 'fail_sound.wav',
    secs: (- 1),
    });
  fail_soundPrac.setVolume(1);
  correction_audio = new sound.Sound({
    win: psychoJS.window,
    value: 'A',
    secs: (- 1),
    });
  correction_audio.setVolume(1);
  // Initialize components for Routine "clericUnlock"
  clericUnlockClock = new util.Clock();
  clericUnlock_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  clericUnlockImg = new visual.ImageStim({
    win : psychoJS.window,
    name : 'clericUnlockImg', units : 'pix', 
    image : '2instructions1.png', mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  clericUnlockMovClock = new util.Clock();
  clericUnlockMov = new visual.MovieStim({
    win: psychoJS.window,
    name: 'clericUnlockMov',
    units: 'pix',
    movie: 'ClericIdleOpening.mp4',
    pos: [0, 0],
    size: [200, 200],
    ori: 0,
    opacity: 1,
    loop: true,
    noAudio: false,
    });
  clericUnlockWav = new sound.Sound({
    win: psychoJS.window,
    value: 'page4.wav',
    secs: (- 1),
    });
  clericUnlockWav.setVolume(1);
  // Initialize components for Routine "InstLoop"
  InstLoopClock = new util.Clock();
  instloop_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  instImg = new visual.ImageStim({
    win : psychoJS.window,
    name : 'instImg', units : 'pix', 
    image : undefined, mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  pageloop = new sound.Sound({
    win: psychoJS.window,
    value: 'A',
    secs: (- 1),
    });
  pageloop.setVolume(1);
  // Initialize components for Routine "blockSetUp"
  blockSetUpClock = new util.Clock();
  score = 0;
  coins = 0;
  blockCount = 0;
  animMov = "Block1.mp4";
  unlockMov = "RogueUnlock1.mp4";
  charUnlockImg = "4instructions1.png";
  charUnlockWav = "rogueUnlock.wav";
  valleyWav = "25valley.wav";
  valleyImg = "5instructions1.png";
  nPlayEnd = 1;
  
  // Initialize components for Routine "LexBlock"
  LexBlockClock = new util.Clock();
  lexCrossStim = new visual.ShapeStim ({
    win: psychoJS.window, name: 'lexCrossStim', 
    vertices: 'cross', size:[0.015, 0.015],
    ori: 0, pos: [0, 0],
    lineWidth: 1, lineColor: new util.Color([1, 1, 1]),
    fillColor: new util.Color([1, 1, 1]),
    opacity: 1, depth: 0, interpolate: true,
  });
  
  word_stim = new visual.TextStim({
    win: psychoJS.window,
    name: 'word_stim',
    text: 'default text',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0], height: 0.045,  wrapWidth: undefined, ori: 0,
    color: new util.Color('white'),  opacity: 1,
    depth: -1.0 
  });
  
  lexCrossStimAfter = new visual.ShapeStim ({
    win: psychoJS.window, name: 'lexCrossStimAfter', 
    vertices: 'cross', size:[0.015, 0.015],
    ori: 0, pos: [0, 0],
    lineWidth: 1, lineColor: new util.Color([1, 1, 1]),
    fillColor: new util.Color([1, 1, 1]),
    opacity: 1, depth: -2, interpolate: true,
  });
  
  lex_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  lexArrowKey = new visual.ImageStim({
    win : psychoJS.window,
    name : 'lexArrowKey', units : undefined, 
    image : 'arrowkey_lex.png', mask : undefined,
    ori : 0, pos : [0, (- 0.2)], size : [0.2, 0.2],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -4.0 
  });
  // Initialize components for Routine "correctSoundCode"
  correctSoundCodeClock = new util.Clock();
  // Initialize components for Routine "correctSound"
  correctSoundClock = new util.Clock();
  coin_sound_2lex = new sound.Sound({
    win: psychoJS.window,
    value: 'coin_sound.wav',
    secs: (- 1),
    });
  coin_sound_2lex.setVolume(0.5);
  lexARrowKeyCorr_2 = new visual.ImageStim({
    win : psychoJS.window,
    name : 'lexARrowKeyCorr_2', units : undefined, 
    image : 'arrowkey_lex.png', mask : undefined,
    ori : 0, pos : [0, (- 0.2)], size : [0.2, 0.2],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  lexCrossStimCSoundLex = new visual.ShapeStim ({
    win: psychoJS.window, name: 'lexCrossStimCSoundLex', 
    vertices: 'cross', size:[0.015, 0.015],
    ori: 0, pos: [0, 0],
    lineWidth: 1, lineColor: new util.Color([1, 1, 1]),
    fillColor: new util.Color([1, 1, 1]),
    opacity: 1, depth: -2, interpolate: true,
  });
  
  // Initialize components for Routine "inCorrSound"
  inCorrSoundClock = new util.Clock();
  coin_sound_3 = new sound.Sound({
    win: psychoJS.window,
    value: 'fail_sound.wav',
    secs: (- 1),
    });
  coin_sound_3.setVolume(1);
  arrowKeyLexCorr = new visual.ImageStim({
    win : psychoJS.window,
    name : 'arrowKeyLexCorr', units : undefined, 
    image : 'arrowkey_lex.png', mask : undefined,
    ori : 0, pos : [0, (- 0.2)], size : [0.2, 0.2],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  lexCrossStim1 = new visual.ShapeStim ({
    win: psychoJS.window, name: 'lexCrossStim1', 
    vertices: 'cross', size:[0.015, 0.015],
    ori: 0, pos: [0, 0],
    lineWidth: 1, lineColor: new util.Color([1, 1, 1]),
    fillColor: new util.Color([1, 1, 1]),
    opacity: 1, depth: -2, interpolate: true,
  });
  
  // Initialize components for Routine "blockAnimLoop"
  blockAnimLoopClock = new util.Clock();
  coin_icon = new visual.ImageStim({
    win : psychoJS.window,
    name : 'coin_icon', units : 'height', 
    image : 'coinicon.png', mask : undefined,
    ori : 0, pos : [(- 0.05), 0.3], size : [0.075, 0.075],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : 0.0 
  });
  correctAnimText = new visual.TextStim({
    win: psychoJS.window,
    name: 'correctAnimText',
    text: '+10!',
    font: 'Arial',
    units: undefined, 
    pos: [0.025, 0.3], height: 0.05,  wrapWidth: undefined, ori: 0,
    color: new util.Color('white'),  opacity: 1,
    depth: -1.0 
  });
  
  animCoinSound = new sound.Sound({
    win: psychoJS.window,
    value: 'coin_sound.wav',
    secs: (- 1),
    });
  animCoinSound.setVolume(0.5);
  // Initialize components for Routine "fixateOnCross"
  fixateOnCrossClock = new util.Clock();
  lexCross_b4 = new visual.ShapeStim ({
    win: psychoJS.window, name: 'lexCross_b4', 
    vertices: 'cross', size:[0.015, 0.015],
    ori: 0, pos: [0, 0],
    lineWidth: 1, lineColor: new util.Color([1, 1, 1]),
    fillColor: new util.Color([1, 1, 1]),
    opacity: 1, depth: 0, interpolate: true,
  });
  
  lexArrowKey_again = new visual.ImageStim({
    win : psychoJS.window,
    name : 'lexArrowKey_again', units : undefined, 
    image : 'arrowkey_lex.png', mask : undefined,
    ori : 0, pos : [0, (- 0.2)], size : [0.2, 0.2],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  // Initialize components for Routine "halfwayBreak"
  halfwayBreakClock = new util.Clock();
  halfway_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  halfwayImg = new visual.ImageStim({
    win : psychoJS.window,
    name : 'halfwayImg', units : 'pix', 
    image : undefined, mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -1.0 
  });
  halfwayWav = new sound.Sound({
    win: psychoJS.window,
    value: 'A',
    secs: (- 1),
    });
  halfwayWav.setVolume(1);
  halfwaySoundWav = new sound.Sound({
    win: psychoJS.window,
    value: 'halfway_sound.wav',
    secs: (- 1),
    });
  halfwaySoundWav.setVolume(1);
  // Initialize components for Routine "InnInst"
  InnInstClock = new util.Clock();
  valleyinst_img = new visual.ImageStim({
    win : psychoJS.window,
    name : 'valleyinst_img', units : 'pix', 
    image : 'inninstructionsS.png', mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : 0.0 
  });
  valley_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  inn_voice = new sound.Sound({
    win: psychoJS.window,
    value: 'inn.wav',
    secs: (- 1),
    });
  inn_voice.setVolume(1);
  // Initialize components for Routine "charUnlock"
  charUnlockClock = new util.Clock();
  CharUn_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  CharUnlockImg = new visual.ImageStim({
    win : psychoJS.window,
    name : 'CharUnlockImg', units : 'pix', 
    image : undefined, mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -2.0 
  });
  amountofGold = new visual.TextStim({
    win: psychoJS.window,
    name: 'amountofGold',
    text: 'default text',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0.175], height: 0.035,  wrapWidth: undefined, ori: 0,
    color: new util.Color('gold'),  opacity: 1,
    depth: -3.0 
  });
  
  charUn_voice = new sound.Sound({
    win: psychoJS.window,
    value: 'A',
    secs: (- 1),
    });
  charUn_voice.setVolume(1);
  // Initialize components for Routine "valleyInst"
  valleyInstClock = new util.Clock();
  valleyinst_image = new visual.ImageStim({
    win : psychoJS.window,
    name : 'valleyinst_image', units : 'pix', 
    image : undefined, mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : 0.0 
  });
  valley_resp1 = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  valley_voice = new sound.Sound({
    win: psychoJS.window,
    value: 'A',
    secs: (- 1),
    });
  valley_voice.setVolume(1);
  // Initialize components for Routine "gameEnd"
  gameEndClock = new util.Clock();
  theEnd_text = new visual.ImageStim({
    win : psychoJS.window,
    name : 'theEnd_text', units : 'pix', 
    image : '6instructions1.png', mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : 0.0 
  });
  endMovClock = new util.Clock();
  endMov = new visual.MovieStim({
    win: psychoJS.window,
    name: 'endMov',
    units: 'pix',
    movie: 'GameEND.mp4',
    pos: [0, 0],
    size: [822, 298],
    ori: 0,
    opacity: 1,
    loop: false,
    noAudio: false,
    });
  EndImg = new visual.ImageStim({
    win : psychoJS.window,
    name : 'EndImg', units : 'pix', 
    image : 'GameENDFinalStaticImage.png', mask : undefined,
    ori : 0, pos : [0, 0], size : [822, 298],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : -2.0 
  });
  end_voice = new sound.Sound({
    win: psychoJS.window,
    value: 'end.wav',
    secs: (- 1),
    });
  end_voice.setVolume(1);
  key_resp_end = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "Credits"
  CreditsClock = new util.Clock();
  Credits_text = new visual.ImageStim({
    win : psychoJS.window,
    name : 'Credits_text', units : 'pix', 
    image : 'Credits.png', mask : undefined,
    ori : 0, pos : [0, 0], size : [1229, 768],
    color : new util.Color([1, 1, 1]), opacity : 1,
    flipHoriz : false, flipVert : false,
    texRes : 128, interpolate : true, depth : 0.0 
  });
  credits_wav = new sound.Sound({
    win: psychoJS.window,
    value: 'Credit.wav',
    secs: (- 1),
    });
  credits_wav.setVolume(1);
  credits_resp = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Create some handy timers
  globalClock = new util.Clock();  // to track the time since experiment started
  routineTimer = new util.CountdownTimer();  // to track time remaining of each (non-slip) routine
  
  return Scheduler.Event.NEXT;
}


var t;
var frameN;
var exSetUpComponents;
function exSetUpRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'exSetUp'-------
    t = 0;
    exSetUpClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    // keep track of which components have finished
    exSetUpComponents = [];
    
    for (const thisComponent of exSetUpComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


var continueRoutine;
function exSetUpRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'exSetUp'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = exSetUpClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of exSetUpComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function exSetUpRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'exSetUp'-------
    for (const thisComponent of exSetUpComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    // the Routine "exSetUp" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var InstLoop1;
var currentLoop;
function InstLoop1LoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  InstLoop1 = new TrialHandler({
    psychoJS: psychoJS,
    nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
    extraInfo: expInfo, originPath: undefined,
    trialList: 'LexPracInst.xlsx',
    seed: undefined, name: 'InstLoop1'
  });
  psychoJS.experiment.addLoop(InstLoop1); // add the loop to the experiment
  currentLoop = InstLoop1;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisInstLoop1 of InstLoop1) {
    const snapshot = InstLoop1.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(InstLoopRoutineBegin(snapshot));
    thisScheduler.add(InstLoopRoutineEachFrame(snapshot));
    thisScheduler.add(InstLoopRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function InstLoop1LoopEnd() {
  psychoJS.experiment.removeLoop(InstLoop1);

  return Scheduler.Event.NEXT;
}


var pracBlock;
function pracBlockLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  pracBlock = new TrialHandler({
    psychoJS: psychoJS,
    nReps: 1, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: 'blockprac.xlsx',
    seed: undefined, name: 'pracBlock'
  });
  psychoJS.experiment.addLoop(pracBlock); // add the loop to the experiment
  currentLoop = pracBlock;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisPracBlock of pracBlock) {
    const snapshot = pracBlock.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(LexBlockRoutineBegin(snapshot));
    thisScheduler.add(LexBlockRoutineEachFrame(snapshot));
    thisScheduler.add(LexBlockRoutineEnd(snapshot));
    thisScheduler.add(pCorrectSoundCodeRoutineBegin(snapshot));
    thisScheduler.add(pCorrectSoundCodeRoutineEachFrame(snapshot));
    thisScheduler.add(pCorrectSoundCodeRoutineEnd(snapshot));
    const pCorrectLoopLoopScheduler = new Scheduler(psychoJS);
    thisScheduler.add(pCorrectLoopLoopBegin, pCorrectLoopLoopScheduler);
    thisScheduler.add(pCorrectLoopLoopScheduler);
    thisScheduler.add(pCorrectLoopLoopEnd);
    const nIncorrectLoopLoopScheduler = new Scheduler(psychoJS);
    thisScheduler.add(nIncorrectLoopLoopBegin, nIncorrectLoopLoopScheduler);
    thisScheduler.add(nIncorrectLoopLoopScheduler);
    thisScheduler.add(nIncorrectLoopLoopEnd);
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


var pCorrectLoop;
function pCorrectLoopLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  pCorrectLoop = new TrialHandler({
    psychoJS: psychoJS,
    nReps: ncorrectSound, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: undefined,
    seed: undefined, name: 'pCorrectLoop'
  });
  psychoJS.experiment.addLoop(pCorrectLoop); // add the loop to the experiment
  currentLoop = pCorrectLoop;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisPCorrectLoop of pCorrectLoop) {
    const snapshot = pCorrectLoop.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(correctSoundPracRoutineBegin(snapshot));
    thisScheduler.add(correctSoundPracRoutineEachFrame(snapshot));
    thisScheduler.add(correctSoundPracRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function pCorrectLoopLoopEnd() {
  psychoJS.experiment.removeLoop(pCorrectLoop);

  return Scheduler.Event.NEXT;
}


var nIncorrectLoop;
function nIncorrectLoopLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  nIncorrectLoop = new TrialHandler({
    psychoJS: psychoJS,
    nReps: nincorrectSound, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: undefined,
    seed: undefined, name: 'nIncorrectLoop'
  });
  psychoJS.experiment.addLoop(nIncorrectLoop); // add the loop to the experiment
  currentLoop = nIncorrectLoop;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisNIncorrectLoop of nIncorrectLoop) {
    const snapshot = nIncorrectLoop.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(incorrectPracRoutineBegin(snapshot));
    thisScheduler.add(incorrectPracRoutineEachFrame(snapshot));
    thisScheduler.add(incorrectPracRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function nIncorrectLoopLoopEnd() {
  psychoJS.experiment.removeLoop(nIncorrectLoop);

  return Scheduler.Event.NEXT;
}


function pracBlockLoopEnd() {
  psychoJS.experiment.removeLoop(pracBlock);

  return Scheduler.Event.NEXT;
}


var InstLoop2;
function InstLoop2LoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  InstLoop2 = new TrialHandler({
    psychoJS: psychoJS,
    nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
    extraInfo: expInfo, originPath: undefined,
    trialList: 'LexPracInst2.xlsx',
    seed: undefined, name: 'InstLoop2'
  });
  psychoJS.experiment.addLoop(InstLoop2); // add the loop to the experiment
  currentLoop = InstLoop2;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisInstLoop2 of InstLoop2) {
    const snapshot = InstLoop2.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(InstLoopRoutineBegin(snapshot));
    thisScheduler.add(InstLoopRoutineEachFrame(snapshot));
    thisScheduler.add(InstLoopRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function InstLoop2LoopEnd() {
  psychoJS.experiment.removeLoop(InstLoop2);

  return Scheduler.Event.NEXT;
}


var threeBlocksLoop;
function threeBlocksLoopLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  threeBlocksLoop = new TrialHandler({
    psychoJS: psychoJS,
    nReps: 1, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: 'randBlock.xlsx',
    seed: undefined, name: 'threeBlocksLoop'
  });
  psychoJS.experiment.addLoop(threeBlocksLoop); // add the loop to the experiment
  currentLoop = threeBlocksLoop;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisThreeBlocksLoop of threeBlocksLoop) {
    const snapshot = threeBlocksLoop.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(blockSetUpRoutineBegin(snapshot));
    thisScheduler.add(blockSetUpRoutineEachFrame(snapshot));
    thisScheduler.add(blockSetUpRoutineEnd(snapshot));
    const blockSelecterLoopScheduler = new Scheduler(psychoJS);
    thisScheduler.add(blockSelecterLoopBegin, blockSelecterLoopScheduler);
    thisScheduler.add(blockSelecterLoopScheduler);
    thisScheduler.add(blockSelecterLoopEnd);
    const charUnlockDecLoopScheduler = new Scheduler(psychoJS);
    thisScheduler.add(charUnlockDecLoopBegin, charUnlockDecLoopScheduler);
    thisScheduler.add(charUnlockDecLoopScheduler);
    thisScheduler.add(charUnlockDecLoopEnd);
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


var blockSelecter;
function blockSelecterLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  blockSelecter = new TrialHandler({
    psychoJS: psychoJS,
    nReps: 1, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: blockNum,
    seed: undefined, name: 'blockSelecter'
  });
  psychoJS.experiment.addLoop(blockSelecter); // add the loop to the experiment
  currentLoop = blockSelecter;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisBlockSelecter of blockSelecter) {
    const snapshot = blockSelecter.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(LexBlockRoutineBegin(snapshot));
    thisScheduler.add(LexBlockRoutineEachFrame(snapshot));
    thisScheduler.add(LexBlockRoutineEnd(snapshot));
    thisScheduler.add(correctSoundCodeRoutineBegin(snapshot));
    thisScheduler.add(correctSoundCodeRoutineEachFrame(snapshot));
    thisScheduler.add(correctSoundCodeRoutineEnd(snapshot));
    const corSoundPlayLoopScheduler = new Scheduler(psychoJS);
    thisScheduler.add(corSoundPlayLoopBegin, corSoundPlayLoopScheduler);
    thisScheduler.add(corSoundPlayLoopScheduler);
    thisScheduler.add(corSoundPlayLoopEnd);
    const incorSoundPlayLoopScheduler = new Scheduler(psychoJS);
    thisScheduler.add(incorSoundPlayLoopBegin, incorSoundPlayLoopScheduler);
    thisScheduler.add(incorSoundPlayLoopScheduler);
    thisScheduler.add(incorSoundPlayLoopEnd);
    const animLoopPlayLoopScheduler = new Scheduler(psychoJS);
    thisScheduler.add(animLoopPlayLoopBegin, animLoopPlayLoopScheduler);
    thisScheduler.add(animLoopPlayLoopScheduler);
    thisScheduler.add(animLoopPlayLoopEnd);
    const halfwayLoopScheduler = new Scheduler(psychoJS);
    thisScheduler.add(halfwayLoopBegin, halfwayLoopScheduler);
    thisScheduler.add(halfwayLoopScheduler);
    thisScheduler.add(halfwayLoopEnd);
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


var corSoundPlay;
function corSoundPlayLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  corSoundPlay = new TrialHandler({
    psychoJS: psychoJS,
    nReps: ncorrectSound, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: undefined,
    seed: undefined, name: 'corSoundPlay'
  });
  psychoJS.experiment.addLoop(corSoundPlay); // add the loop to the experiment
  currentLoop = corSoundPlay;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisCorSoundPlay of corSoundPlay) {
    const snapshot = corSoundPlay.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(correctSoundRoutineBegin(snapshot));
    thisScheduler.add(correctSoundRoutineEachFrame(snapshot));
    thisScheduler.add(correctSoundRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function corSoundPlayLoopEnd() {
  psychoJS.experiment.removeLoop(corSoundPlay);

  return Scheduler.Event.NEXT;
}


var incorSoundPlay;
function incorSoundPlayLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  incorSoundPlay = new TrialHandler({
    psychoJS: psychoJS,
    nReps: nincorrectSound, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: undefined,
    seed: undefined, name: 'incorSoundPlay'
  });
  psychoJS.experiment.addLoop(incorSoundPlay); // add the loop to the experiment
  currentLoop = incorSoundPlay;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisIncorSoundPlay of incorSoundPlay) {
    const snapshot = incorSoundPlay.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(inCorrSoundRoutineBegin(snapshot));
    thisScheduler.add(inCorrSoundRoutineEachFrame(snapshot));
    thisScheduler.add(inCorrSoundRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function incorSoundPlayLoopEnd() {
  psychoJS.experiment.removeLoop(incorSoundPlay);

  return Scheduler.Event.NEXT;
}


var animLoopPlay;
function animLoopPlayLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  animLoopPlay = new TrialHandler({
    psychoJS: psychoJS,
    nReps: ncorrectAnim, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: undefined,
    seed: undefined, name: 'animLoopPlay'
  });
  psychoJS.experiment.addLoop(animLoopPlay); // add the loop to the experiment
  currentLoop = animLoopPlay;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisAnimLoopPlay of animLoopPlay) {
    const snapshot = animLoopPlay.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(blockAnimLoopRoutineBegin(snapshot));
    thisScheduler.add(blockAnimLoopRoutineEachFrame(snapshot));
    thisScheduler.add(blockAnimLoopRoutineEnd(snapshot));
    thisScheduler.add(fixateOnCrossRoutineBegin(snapshot));
    thisScheduler.add(fixateOnCrossRoutineEachFrame(snapshot));
    thisScheduler.add(fixateOnCrossRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function animLoopPlayLoopEnd() {
  psychoJS.experiment.removeLoop(animLoopPlay);

  return Scheduler.Event.NEXT;
}


var halfway;
function halfwayLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  halfway = new TrialHandler({
    psychoJS: psychoJS,
    nReps: nhalfway, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: undefined,
    seed: undefined, name: 'halfway'
  });
  psychoJS.experiment.addLoop(halfway); // add the loop to the experiment
  currentLoop = halfway;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisHalfway of halfway) {
    const snapshot = halfway.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(halfwayBreakRoutineBegin(snapshot));
    thisScheduler.add(halfwayBreakRoutineEachFrame(snapshot));
    thisScheduler.add(halfwayBreakRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function halfwayLoopEnd() {
  psychoJS.experiment.removeLoop(halfway);

  return Scheduler.Event.NEXT;
}


function blockSelecterLoopEnd() {
  psychoJS.experiment.removeLoop(blockSelecter);

  return Scheduler.Event.NEXT;
}


var charUnlockDec;
function charUnlockDecLoopBegin(thisScheduler) {
  // set up handler to look after randomisation of conditions etc
  charUnlockDec = new TrialHandler({
    psychoJS: psychoJS,
    nReps: nPlayEnd, method: TrialHandler.Method.RANDOM,
    extraInfo: expInfo, originPath: undefined,
    trialList: undefined,
    seed: undefined, name: 'charUnlockDec'
  });
  psychoJS.experiment.addLoop(charUnlockDec); // add the loop to the experiment
  currentLoop = charUnlockDec;  // we're now the current loop

  // Schedule all the trials in the trialList:
  for (const thisCharUnlockDec of charUnlockDec) {
    const snapshot = charUnlockDec.getSnapshot();
    thisScheduler.add(importConditions(snapshot));
    thisScheduler.add(InnInstRoutineBegin(snapshot));
    thisScheduler.add(InnInstRoutineEachFrame(snapshot));
    thisScheduler.add(InnInstRoutineEnd(snapshot));
    thisScheduler.add(charUnlockRoutineBegin(snapshot));
    thisScheduler.add(charUnlockRoutineEachFrame(snapshot));
    thisScheduler.add(charUnlockRoutineEnd(snapshot));
    thisScheduler.add(valleyInstRoutineBegin(snapshot));
    thisScheduler.add(valleyInstRoutineEachFrame(snapshot));
    thisScheduler.add(valleyInstRoutineEnd(snapshot));
    thisScheduler.add(endLoopIteration(thisScheduler, snapshot));
  }

  return Scheduler.Event.NEXT;
}


function charUnlockDecLoopEnd() {
  psychoJS.experiment.removeLoop(charUnlockDec);

  return Scheduler.Event.NEXT;
}


function threeBlocksLoopLoopEnd() {
  psychoJS.experiment.removeLoop(threeBlocksLoop);

  return Scheduler.Event.NEXT;
}


var _instloop_resp_allKeys;
var InstLoopComponents;
function InstLoopRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'InstLoop'-------
    t = 0;
    InstLoopClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    instloop_resp.keys = undefined;
    instloop_resp.rt = undefined;
    _instloop_resp_allKeys = [];
    instImg.setImage(instructionpng);
    pageloop = new sound.Sound({
    win: psychoJS.window,
    value: soundwav,
    secs: -1,
    });
    pageloop.setVolume(1);
    // keep track of which components have finished
    InstLoopComponents = [];
    InstLoopComponents.push(instloop_resp);
    InstLoopComponents.push(instImg);
    InstLoopComponents.push(pageloop);
    
    for (const thisComponent of InstLoopComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function InstLoopRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'InstLoop'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = InstLoopClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *instloop_resp* updates
    if (t >= 1 && instloop_resp.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      instloop_resp.tStart = t;  // (not accounting for frame time here)
      instloop_resp.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { instloop_resp.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { instloop_resp.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { instloop_resp.clearEvents(); });
    }

    if (instloop_resp.status === PsychoJS.Status.STARTED) {
      let theseKeys = instloop_resp.getKeys({keyList: ['space'], waitRelease: false});
      _instloop_resp_allKeys = _instloop_resp_allKeys.concat(theseKeys);
      if (_instloop_resp_allKeys.length > 0) {
        instloop_resp.keys = _instloop_resp_allKeys[_instloop_resp_allKeys.length - 1].name;  // just the last key pressed
        instloop_resp.rt = _instloop_resp_allKeys[_instloop_resp_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *instImg* updates
    if (t >= 0.0 && instImg.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      instImg.tStart = t;  // (not accounting for frame time here)
      instImg.frameNStart = frameN;  // exact frame index
      
      instImg.setAutoDraw(true);
    }

    // start/stop pageloop
    if (t >= 0.3 && pageloop.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pageloop.tStart = t;  // (not accounting for frame time here)
      pageloop.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ pageloop.play(); });  // screen flip
      pageloop.status = PsychoJS.Status.STARTED;
    }
    if (t >= (pageloop.getDuration() + pageloop.tStart)     && pageloop.status === PsychoJS.Status.STARTED) {
      pageloop.stop();  // stop the sound (if longer than duration)
      pageloop.status = PsychoJS.Status.FINISHED;
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of InstLoopComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function InstLoopRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'InstLoop'-------
    for (const thisComponent of InstLoopComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('instloop_resp.keys', instloop_resp.keys);
    if (typeof instloop_resp.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('instloop_resp.rt', instloop_resp.rt);
        routineTimer.reset();
        }
    
    instloop_resp.stop();
    pageloop.stop();  // ensure sound has stopped at end of routine
    // the Routine "InstLoop" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _lex_resp_allKeys;
var LexBlockComponents;
function LexBlockRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'LexBlock'-------
    t = 0;
    LexBlockClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    word_stim.setText(words);
    lex_resp.keys = undefined;
    lex_resp.rt = undefined;
    _lex_resp_allKeys = [];
    // keep track of which components have finished
    LexBlockComponents = [];
    LexBlockComponents.push(lexCrossStim);
    LexBlockComponents.push(word_stim);
    LexBlockComponents.push(lexCrossStimAfter);
    LexBlockComponents.push(lex_resp);
    LexBlockComponents.push(lexArrowKey);
    
    for (const thisComponent of LexBlockComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


var frameRemains;
function LexBlockRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'LexBlock'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = LexBlockClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *lexCrossStim* updates
    if (t >= 0.0 && lexCrossStim.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexCrossStim.tStart = t;  // (not accounting for frame time here)
      lexCrossStim.frameNStart = frameN;  // exact frame index
      
      lexCrossStim.setAutoDraw(true);
    }

    frameRemains = 0.0 + 0.5 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (lexCrossStim.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      lexCrossStim.setAutoDraw(false);
    }
    
    // *word_stim* updates
    if (t >= 0.5 && word_stim.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      word_stim.tStart = t;  // (not accounting for frame time here)
      word_stim.frameNStart = frameN;  // exact frame index
      
      word_stim.setAutoDraw(true);
    }

    frameRemains = 0.5 + 0.35 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (word_stim.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      word_stim.setAutoDraw(false);
    }
    
    // *lexCrossStimAfter* updates
    if (t >= 0.85 && lexCrossStimAfter.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexCrossStimAfter.tStart = t;  // (not accounting for frame time here)
      lexCrossStimAfter.frameNStart = frameN;  // exact frame index
      
      lexCrossStimAfter.setAutoDraw(true);
    }

    
    // *lex_resp* updates
    if (t >= 0.5 && lex_resp.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lex_resp.tStart = t;  // (not accounting for frame time here)
      lex_resp.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { lex_resp.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { lex_resp.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { lex_resp.clearEvents(); });
    }

    if (lex_resp.status === PsychoJS.Status.STARTED) {
      let theseKeys = lex_resp.getKeys({keyList: ['left', 'right'], waitRelease: false});
      _lex_resp_allKeys = _lex_resp_allKeys.concat(theseKeys);
      if (_lex_resp_allKeys.length > 0) {
        lex_resp.keys = _lex_resp_allKeys[_lex_resp_allKeys.length - 1].name;  // just the last key pressed
        lex_resp.rt = _lex_resp_allKeys[_lex_resp_allKeys.length - 1].rt;
        // was this correct?
        if (lex_resp.keys == correctAns) {
            lex_resp.corr = 1;
        } else {
            lex_resp.corr = 0;
        }
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *lexArrowKey* updates
    if (t >= 0.0 && lexArrowKey.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexArrowKey.tStart = t;  // (not accounting for frame time here)
      lexArrowKey.frameNStart = frameN;  // exact frame index
      
      lexArrowKey.setAutoDraw(true);
    }

    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of LexBlockComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function LexBlockRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'LexBlock'-------
    for (const thisComponent of LexBlockComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    // was no response the correct answer?!
    if (lex_resp.keys === undefined) {
      if (['None','none',undefined].includes(correctAns)) {
         lex_resp.corr = 1;  // correct non-response
      } else {
         lex_resp.corr = 0;  // failed to respond (incorrectly)
      }
    }
    // store data for thisExp (ExperimentHandler)
    psychoJS.experiment.addData('lex_resp.keys', lex_resp.keys);
    psychoJS.experiment.addData('lex_resp.corr', lex_resp.corr);
    if (typeof lex_resp.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('lex_resp.rt', lex_resp.rt);
        routineTimer.reset();
        }
    
    lex_resp.stop();
    // the Routine "LexBlock" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var pCorrectSoundCodeComponents;
function pCorrectSoundCodeRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'pCorrectSoundCode'-------
    t = 0;
    pCorrectSoundCodeClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    if ((lex_resp.corr === 1)) {
        ncorrectSound = 1;
        nincorrectSound = 0;
    } else {
        ncorrectSound = 0;
        nincorrectSound = 1;
    }
    
    // keep track of which components have finished
    pCorrectSoundCodeComponents = [];
    
    for (const thisComponent of pCorrectSoundCodeComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function pCorrectSoundCodeRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'pCorrectSoundCode'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = pCorrectSoundCodeClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of pCorrectSoundCodeComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function pCorrectSoundCodeRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'pCorrectSoundCode'-------
    for (const thisComponent of pCorrectSoundCodeComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    // the Routine "pCorrectSoundCode" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var correctSoundPracComponents;
function correctSoundPracRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'correctSoundPrac'-------
    t = 0;
    correctSoundPracClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    coin_sound_2.setVolume(0.5);
    // keep track of which components have finished
    correctSoundPracComponents = [];
    correctSoundPracComponents.push(coin_sound_2);
    correctSoundPracComponents.push(triArrowKeyCorr_2);
    correctSoundPracComponents.push(lexCrossStimCSound);
    
    for (const thisComponent of correctSoundPracComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function correctSoundPracRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'correctSoundPrac'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = correctSoundPracClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // start/stop coin_sound_2
    if (t >= 0.0 && coin_sound_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      coin_sound_2.tStart = t;  // (not accounting for frame time here)
      coin_sound_2.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ coin_sound_2.play(); });  // screen flip
      coin_sound_2.status = PsychoJS.Status.STARTED;
    }
    if (t >= (coin_sound_2.getDuration() + coin_sound_2.tStart)     && coin_sound_2.status === PsychoJS.Status.STARTED) {
      coin_sound_2.stop();  // stop the sound (if longer than duration)
      coin_sound_2.status = PsychoJS.Status.FINISHED;
    }
    
    // *triArrowKeyCorr_2* updates
    if (t >= 0.0 && triArrowKeyCorr_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      triArrowKeyCorr_2.tStart = t;  // (not accounting for frame time here)
      triArrowKeyCorr_2.frameNStart = frameN;  // exact frame index
      
      triArrowKeyCorr_2.setAutoDraw(true);
    }

    frameRemains = 0.0 + 0.5 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (triArrowKeyCorr_2.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      triArrowKeyCorr_2.setAutoDraw(false);
    }
    
    // *lexCrossStimCSound* updates
    if (t >= 0.0 && lexCrossStimCSound.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexCrossStimCSound.tStart = t;  // (not accounting for frame time here)
      lexCrossStimCSound.frameNStart = frameN;  // exact frame index
      
      lexCrossStimCSound.setAutoDraw(true);
    }

    frameRemains = 0.0 + 0.5 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (lexCrossStimCSound.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      lexCrossStimCSound.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of correctSoundPracComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function correctSoundPracRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'correctSoundPrac'-------
    for (const thisComponent of correctSoundPracComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    coin_sound_2.stop();  // ensure sound has stopped at end of routine
    // the Routine "correctSoundPrac" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _incorrectPrac_resp_allKeys;
var incorrectPracComponents;
function incorrectPracRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'incorrectPrac'-------
    t = 0;
    incorrectPracClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    incorrectPrac_resp.keys = undefined;
    incorrectPrac_resp.rt = undefined;
    _incorrectPrac_resp_allKeys = [];
    incorrect_image.setImage(corrpng);
    fail_soundPrac.setVolume(1);
    correction_audio = new sound.Sound({
    win: psychoJS.window,
    value: corrwav,
    secs: -1,
    });
    correction_audio.setVolume(1);
    // keep track of which components have finished
    incorrectPracComponents = [];
    incorrectPracComponents.push(incorrectPrac_resp);
    incorrectPracComponents.push(incorrect_image);
    incorrectPracComponents.push(fail_soundPrac);
    incorrectPracComponents.push(correction_audio);
    
    for (const thisComponent of incorrectPracComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function incorrectPracRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'incorrectPrac'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = incorrectPracClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *incorrectPrac_resp* updates
    if (t >= 0.0 && incorrectPrac_resp.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      incorrectPrac_resp.tStart = t;  // (not accounting for frame time here)
      incorrectPrac_resp.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { incorrectPrac_resp.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { incorrectPrac_resp.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { incorrectPrac_resp.clearEvents(); });
    }

    if (incorrectPrac_resp.status === PsychoJS.Status.STARTED) {
      let theseKeys = incorrectPrac_resp.getKeys({keyList: ['space'], waitRelease: false});
      _incorrectPrac_resp_allKeys = _incorrectPrac_resp_allKeys.concat(theseKeys);
      if (_incorrectPrac_resp_allKeys.length > 0) {
        incorrectPrac_resp.keys = _incorrectPrac_resp_allKeys[_incorrectPrac_resp_allKeys.length - 1].name;  // just the last key pressed
        incorrectPrac_resp.rt = _incorrectPrac_resp_allKeys[_incorrectPrac_resp_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *incorrect_image* updates
    if (t >= 0.0 && incorrect_image.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      incorrect_image.tStart = t;  // (not accounting for frame time here)
      incorrect_image.frameNStart = frameN;  // exact frame index
      
      incorrect_image.setAutoDraw(true);
    }

    // start/stop fail_soundPrac
    if (t >= 0.0 && fail_soundPrac.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      fail_soundPrac.tStart = t;  // (not accounting for frame time here)
      fail_soundPrac.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ fail_soundPrac.play(); });  // screen flip
      fail_soundPrac.status = PsychoJS.Status.STARTED;
    }
    if (t >= (fail_soundPrac.getDuration() + fail_soundPrac.tStart)     && fail_soundPrac.status === PsychoJS.Status.STARTED) {
      fail_soundPrac.stop();  // stop the sound (if longer than duration)
      fail_soundPrac.status = PsychoJS.Status.FINISHED;
    }
    // start/stop correction_audio
    if (t >= 0.0 && correction_audio.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      correction_audio.tStart = t;  // (not accounting for frame time here)
      correction_audio.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ correction_audio.play(); });  // screen flip
      correction_audio.status = PsychoJS.Status.STARTED;
    }
    if (t >= (correction_audio.getDuration() + correction_audio.tStart)     && correction_audio.status === PsychoJS.Status.STARTED) {
      correction_audio.stop();  // stop the sound (if longer than duration)
      correction_audio.status = PsychoJS.Status.FINISHED;
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of incorrectPracComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function incorrectPracRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'incorrectPrac'-------
    for (const thisComponent of incorrectPracComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('incorrectPrac_resp.keys', incorrectPrac_resp.keys);
    if (typeof incorrectPrac_resp.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('incorrectPrac_resp.rt', incorrectPrac_resp.rt);
        routineTimer.reset();
        }
    
    incorrectPrac_resp.stop();
    fail_soundPrac.stop();  // ensure sound has stopped at end of routine
    correction_audio.stop();  // ensure sound has stopped at end of routine
    // the Routine "incorrectPrac" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _clericUnlock_resp_allKeys;
var clericUnlockComponents;
function clericUnlockRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'clericUnlock'-------
    t = 0;
    clericUnlockClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    clericUnlock_resp.keys = undefined;
    clericUnlock_resp.rt = undefined;
    _clericUnlock_resp_allKeys = [];
    clericUnlockWav.setVolume(1);
    // keep track of which components have finished
    clericUnlockComponents = [];
    clericUnlockComponents.push(clericUnlock_resp);
    clericUnlockComponents.push(clericUnlockImg);
    clericUnlockComponents.push(clericUnlockMov);
    clericUnlockComponents.push(clericUnlockWav);
    
    for (const thisComponent of clericUnlockComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function clericUnlockRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'clericUnlock'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = clericUnlockClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *clericUnlock_resp* updates
    if (t >= 0.0 && clericUnlock_resp.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      clericUnlock_resp.tStart = t;  // (not accounting for frame time here)
      clericUnlock_resp.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { clericUnlock_resp.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { clericUnlock_resp.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { clericUnlock_resp.clearEvents(); });
    }

    if (clericUnlock_resp.status === PsychoJS.Status.STARTED) {
      let theseKeys = clericUnlock_resp.getKeys({keyList: ['space'], waitRelease: false});
      _clericUnlock_resp_allKeys = _clericUnlock_resp_allKeys.concat(theseKeys);
      if (_clericUnlock_resp_allKeys.length > 0) {
        clericUnlock_resp.keys = _clericUnlock_resp_allKeys[_clericUnlock_resp_allKeys.length - 1].name;  // just the last key pressed
        clericUnlock_resp.rt = _clericUnlock_resp_allKeys[_clericUnlock_resp_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *clericUnlockImg* updates
    if (t >= 0.0 && clericUnlockImg.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      clericUnlockImg.tStart = t;  // (not accounting for frame time here)
      clericUnlockImg.frameNStart = frameN;  // exact frame index
      
      clericUnlockImg.setAutoDraw(true);
    }

    
    // *clericUnlockMov* updates
    if (t >= 0.0 && clericUnlockMov.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      clericUnlockMov.tStart = t;  // (not accounting for frame time here)
      clericUnlockMov.frameNStart = frameN;  // exact frame index
      
      clericUnlockMov.setAutoDraw(true);
    }

    // start/stop clericUnlockWav
    if (t >= 0.0 && clericUnlockWav.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      clericUnlockWav.tStart = t;  // (not accounting for frame time here)
      clericUnlockWav.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ clericUnlockWav.play(); });  // screen flip
      clericUnlockWav.status = PsychoJS.Status.STARTED;
    }
    if (t >= (clericUnlockWav.getDuration() + clericUnlockWav.tStart)     && clericUnlockWav.status === PsychoJS.Status.STARTED) {
      clericUnlockWav.stop();  // stop the sound (if longer than duration)
      clericUnlockWav.status = PsychoJS.Status.FINISHED;
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of clericUnlockComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function clericUnlockRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'clericUnlock'-------
    for (const thisComponent of clericUnlockComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('clericUnlock_resp.keys', clericUnlock_resp.keys);
    if (typeof clericUnlock_resp.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('clericUnlock_resp.rt', clericUnlock_resp.rt);
        routineTimer.reset();
        }
    
    clericUnlock_resp.stop();
    clericUnlockWav.stop();  // ensure sound has stopped at end of routine
    // the Routine "clericUnlock" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var trialCount;
var blockSetUpComponents;
function blockSetUpRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'blockSetUp'-------
    t = 0;
    blockSetUpClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    trialCount = 0;
    blockCount = (blockCount + 1);
    if ((blockCount === 2)) {
        animMov = "Block2.mp4";
        unlockMov = "KnightUnlock2.mp4";
        charUnlockImg = "4instructions2.png";
        charUnlockWav = "knightUnlock.wav";
        valleyWav = "34valley.wav";
        valleyImg = "5instructions2.png";
    } else {
        if ((blockCount === 3)) {
            animMov = "Block3.mp4";
            nPlayEnd = 0;
            valleyWav = "34valley.wav";
            valleyImg = "5instructions2.png";
        }
    }
    
    // keep track of which components have finished
    blockSetUpComponents = [];
    
    for (const thisComponent of blockSetUpComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function blockSetUpRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'blockSetUp'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = blockSetUpClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of blockSetUpComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function blockSetUpRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'blockSetUp'-------
    for (const thisComponent of blockSetUpComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    // the Routine "blockSetUp" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var ncorrectAnim;
var correctSoundCodeComponents;
function correctSoundCodeRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'correctSoundCode'-------
    t = 0;
    correctSoundCodeClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    trialCount = (trialCount + 1);
    ncorrectSound = 0;
    nincorrectSound = 0;
    ncorrectAnim = 0;
    nhalfway = 0;
    if ((lex_resp.corr === 1)) {
        ncorrectSound = 1;
        score = (score + 1);
        if (((score % 10) === 0)) {
            ncorrectSound = 0;
            ncorrectAnim = 1;
            coins = (coins + 10);
        }
        if ((trialCount === 42)) {
            ncorrectSound = 0;
            ncorrectAnim = 0;
            nhalfway = 1;
        }
        if ((trialCount === 84)) {
            ncorrectSound = 0;
            ncorrectAnim = 0;
            trialCount = 0;
        }
    } else {
        nincorrectSound = 1;
        if ((trialCount === 42)) {
            nincorrectSound = 0;
            nhalfway = 1;
        }
        if ((trialCount === 84)) {
            nincorrectSound = 0;
            trialCount = 0;
        }
    }
    
    // keep track of which components have finished
    correctSoundCodeComponents = [];
    
    for (const thisComponent of correctSoundCodeComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function correctSoundCodeRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'correctSoundCode'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = correctSoundCodeClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of correctSoundCodeComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function correctSoundCodeRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'correctSoundCode'-------
    for (const thisComponent of correctSoundCodeComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    // the Routine "correctSoundCode" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var correctSoundComponents;
function correctSoundRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'correctSound'-------
    t = 0;
    correctSoundClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    coin_sound_2lex.setVolume(0.5);
    // keep track of which components have finished
    correctSoundComponents = [];
    correctSoundComponents.push(coin_sound_2lex);
    correctSoundComponents.push(lexARrowKeyCorr_2);
    correctSoundComponents.push(lexCrossStimCSoundLex);
    
    for (const thisComponent of correctSoundComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function correctSoundRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'correctSound'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = correctSoundClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // start/stop coin_sound_2lex
    if (t >= 0.0 && coin_sound_2lex.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      coin_sound_2lex.tStart = t;  // (not accounting for frame time here)
      coin_sound_2lex.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ coin_sound_2lex.play(); });  // screen flip
      coin_sound_2lex.status = PsychoJS.Status.STARTED;
    }
    if (t >= (coin_sound_2lex.getDuration() + coin_sound_2lex.tStart)     && coin_sound_2lex.status === PsychoJS.Status.STARTED) {
      coin_sound_2lex.stop();  // stop the sound (if longer than duration)
      coin_sound_2lex.status = PsychoJS.Status.FINISHED;
    }
    
    // *lexARrowKeyCorr_2* updates
    if (t >= 0.0 && lexARrowKeyCorr_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexARrowKeyCorr_2.tStart = t;  // (not accounting for frame time here)
      lexARrowKeyCorr_2.frameNStart = frameN;  // exact frame index
      
      lexARrowKeyCorr_2.setAutoDraw(true);
    }

    frameRemains = 0.0 + 0.5 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (lexARrowKeyCorr_2.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      lexARrowKeyCorr_2.setAutoDraw(false);
    }
    
    // *lexCrossStimCSoundLex* updates
    if (t >= 0.0 && lexCrossStimCSoundLex.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexCrossStimCSoundLex.tStart = t;  // (not accounting for frame time here)
      lexCrossStimCSoundLex.frameNStart = frameN;  // exact frame index
      
      lexCrossStimCSoundLex.setAutoDraw(true);
    }

    frameRemains = 0.0 + 0.5 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (lexCrossStimCSoundLex.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      lexCrossStimCSoundLex.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of correctSoundComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function correctSoundRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'correctSound'-------
    for (const thisComponent of correctSoundComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    coin_sound_2lex.stop();  // ensure sound has stopped at end of routine
    // the Routine "correctSound" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var inCorrSoundComponents;
function inCorrSoundRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'inCorrSound'-------
    t = 0;
    inCorrSoundClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    coin_sound_3.setVolume(1);
    // keep track of which components have finished
    inCorrSoundComponents = [];
    inCorrSoundComponents.push(coin_sound_3);
    inCorrSoundComponents.push(arrowKeyLexCorr);
    inCorrSoundComponents.push(lexCrossStim1);
    
    for (const thisComponent of inCorrSoundComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function inCorrSoundRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'inCorrSound'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = inCorrSoundClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // start/stop coin_sound_3
    if (t >= 0.0 && coin_sound_3.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      coin_sound_3.tStart = t;  // (not accounting for frame time here)
      coin_sound_3.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ coin_sound_3.play(); });  // screen flip
      coin_sound_3.status = PsychoJS.Status.STARTED;
    }
    if (t >= (coin_sound_3.getDuration() + coin_sound_3.tStart)     && coin_sound_3.status === PsychoJS.Status.STARTED) {
      coin_sound_3.stop();  // stop the sound (if longer than duration)
      coin_sound_3.status = PsychoJS.Status.FINISHED;
    }
    
    // *arrowKeyLexCorr* updates
    if (t >= 0.0 && arrowKeyLexCorr.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      arrowKeyLexCorr.tStart = t;  // (not accounting for frame time here)
      arrowKeyLexCorr.frameNStart = frameN;  // exact frame index
      
      arrowKeyLexCorr.setAutoDraw(true);
    }

    frameRemains = 0.0 + 0.5 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (arrowKeyLexCorr.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      arrowKeyLexCorr.setAutoDraw(false);
    }
    
    // *lexCrossStim1* updates
    if (t >= 0.0 && lexCrossStim1.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexCrossStim1.tStart = t;  // (not accounting for frame time here)
      lexCrossStim1.frameNStart = frameN;  // exact frame index
      
      lexCrossStim1.setAutoDraw(true);
    }

    frameRemains = 0.0 + 0.5 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (lexCrossStim1.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      lexCrossStim1.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of inCorrSoundComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function inCorrSoundRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'inCorrSound'-------
    for (const thisComponent of inCorrSoundComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    coin_sound_3.stop();  // ensure sound has stopped at end of routine
    // the Routine "inCorrSound" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var correctAnimMovClock;
var correctAnimMov;
var blockAnimLoopComponents;
function blockAnimLoopRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'blockAnimLoop'-------
    t = 0;
    blockAnimLoopClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    animCoinSound.setVolume(0.5);
    correctAnimMovClock = new util.Clock();
    correctAnimMov = new visual.MovieStim({
      win: psychoJS.window,
      name: 'correctAnimMov',
      units: 'height',
      movie: animMov,
      pos: [0, 0],
      size: [1.3152, 0.4776],
      ori: 0,
      opacity: 1,
      loop: false,
      noAudio: false,
      });
    correctAnimMov.play();
    
    // keep track of which components have finished
    blockAnimLoopComponents = [];
    blockAnimLoopComponents.push(coin_icon);
    blockAnimLoopComponents.push(correctAnimText);
    blockAnimLoopComponents.push(animCoinSound);
    blockAnimLoopComponents.push(correctAnimMov);
    
    for (const thisComponent of blockAnimLoopComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function blockAnimLoopRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'blockAnimLoop'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = blockAnimLoopClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *coin_icon* updates
    if (t >= 0.0 && coin_icon.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      coin_icon.tStart = t;  // (not accounting for frame time here)
      coin_icon.frameNStart = frameN;  // exact frame index
      
      coin_icon.setAutoDraw(true);
    }

    frameRemains = 0.0 + 1.7 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (coin_icon.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      coin_icon.setAutoDraw(false);
    }
    
    // *correctAnimText* updates
    if (t >= 0.0 && correctAnimText.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      correctAnimText.tStart = t;  // (not accounting for frame time here)
      correctAnimText.frameNStart = frameN;  // exact frame index
      
      correctAnimText.setAutoDraw(true);
    }

    frameRemains = 0.0 + 1.7 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (correctAnimText.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      correctAnimText.setAutoDraw(false);
    }
    // start/stop animCoinSound
    if (t >= 0.0 && animCoinSound.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      animCoinSound.tStart = t;  // (not accounting for frame time here)
      animCoinSound.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ animCoinSound.play(); });  // screen flip
      animCoinSound.status = PsychoJS.Status.STARTED;
    }
    if (t >= (animCoinSound.getDuration() + animCoinSound.tStart)     && animCoinSound.status === PsychoJS.Status.STARTED) {
      animCoinSound.stop();  // stop the sound (if longer than duration)
      animCoinSound.status = PsychoJS.Status.FINISHED;
    }
    
    // *correctAnimMov* updates
    if (t >= 0.0 && correctAnimMov.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      correctAnimMov.tStart = t;  // (not accounting for frame time here)
      correctAnimMov.frameNStart = frameN;  // exact frame index
      
      correctAnimMov.setAutoDraw(true);
    }

    frameRemains = 0.0 + undefined - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (correctAnimMov.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      correctAnimMov.setAutoDraw(false);
    }

    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of blockAnimLoopComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function blockAnimLoopRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'blockAnimLoop'-------
    for (const thisComponent of blockAnimLoopComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    animCoinSound.stop();  // ensure sound has stopped at end of routine
    correctAnimMov.stop();
    
    // the Routine "blockAnimLoop" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var fixateOnCrossComponents;
function fixateOnCrossRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'fixateOnCross'-------
    t = 0;
    fixateOnCrossClock.reset(); // clock
    frameN = -1;
    routineTimer.add(1.000000);
    // update component parameters for each repeat
    // keep track of which components have finished
    fixateOnCrossComponents = [];
    fixateOnCrossComponents.push(lexCross_b4);
    fixateOnCrossComponents.push(lexArrowKey_again);
    
    for (const thisComponent of fixateOnCrossComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function fixateOnCrossRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'fixateOnCross'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = fixateOnCrossClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *lexCross_b4* updates
    if (t >= 0.0 && lexCross_b4.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexCross_b4.tStart = t;  // (not accounting for frame time here)
      lexCross_b4.frameNStart = frameN;  // exact frame index
      
      lexCross_b4.setAutoDraw(true);
    }

    frameRemains = 0.0 + 1.0 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (lexCross_b4.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      lexCross_b4.setAutoDraw(false);
    }
    
    // *lexArrowKey_again* updates
    if (t >= 0.0 && lexArrowKey_again.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      lexArrowKey_again.tStart = t;  // (not accounting for frame time here)
      lexArrowKey_again.frameNStart = frameN;  // exact frame index
      
      lexArrowKey_again.setAutoDraw(true);
    }

    frameRemains = 0.0 + 1.0 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (lexArrowKey_again.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      lexArrowKey_again.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of fixateOnCrossComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine && routineTimer.getTime() > 0) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function fixateOnCrossRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'fixateOnCross'-------
    for (const thisComponent of fixateOnCrossComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    return Scheduler.Event.NEXT;
  };
}


var _halfway_resp_allKeys;
var halfwayBreakComponents;
function halfwayBreakRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'halfwayBreak'-------
    t = 0;
    halfwayBreakClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    halfway_resp.keys = undefined;
    halfway_resp.rt = undefined;
    _halfway_resp_allKeys = [];
    halfwayImg.setImage('8instructions1.png');
    halfwayWav = new sound.Sound({
    win: psychoJS.window,
    value: 'halfway.wav',
    secs: -1,
    });
    halfwayWav.setVolume(1);
    halfwaySoundWav.setVolume(1);
    // keep track of which components have finished
    halfwayBreakComponents = [];
    halfwayBreakComponents.push(halfway_resp);
    halfwayBreakComponents.push(halfwayImg);
    halfwayBreakComponents.push(halfwayWav);
    halfwayBreakComponents.push(halfwaySoundWav);
    
    for (const thisComponent of halfwayBreakComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function halfwayBreakRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'halfwayBreak'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = halfwayBreakClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *halfway_resp* updates
    if (t >= 1 && halfway_resp.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      halfway_resp.tStart = t;  // (not accounting for frame time here)
      halfway_resp.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { halfway_resp.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { halfway_resp.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { halfway_resp.clearEvents(); });
    }

    if (halfway_resp.status === PsychoJS.Status.STARTED) {
      let theseKeys = halfway_resp.getKeys({keyList: ['space'], waitRelease: false});
      _halfway_resp_allKeys = _halfway_resp_allKeys.concat(theseKeys);
      if (_halfway_resp_allKeys.length > 0) {
        halfway_resp.keys = _halfway_resp_allKeys[_halfway_resp_allKeys.length - 1].name;  // just the last key pressed
        halfway_resp.rt = _halfway_resp_allKeys[_halfway_resp_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *halfwayImg* updates
    if (t >= 0.0 && halfwayImg.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      halfwayImg.tStart = t;  // (not accounting for frame time here)
      halfwayImg.frameNStart = frameN;  // exact frame index
      
      halfwayImg.setAutoDraw(true);
    }

    // start/stop halfwayWav
    if (t >= 0.3 && halfwayWav.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      halfwayWav.tStart = t;  // (not accounting for frame time here)
      halfwayWav.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ halfwayWav.play(); });  // screen flip
      halfwayWav.status = PsychoJS.Status.STARTED;
    }
    if (t >= (halfwayWav.getDuration() + halfwayWav.tStart)     && halfwayWav.status === PsychoJS.Status.STARTED) {
      halfwayWav.stop();  // stop the sound (if longer than duration)
      halfwayWav.status = PsychoJS.Status.FINISHED;
    }
    // start/stop halfwaySoundWav
    if (t >= 0.0 && halfwaySoundWav.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      halfwaySoundWav.tStart = t;  // (not accounting for frame time here)
      halfwaySoundWav.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ halfwaySoundWav.play(); });  // screen flip
      halfwaySoundWav.status = PsychoJS.Status.STARTED;
    }
    if (t >= (halfwaySoundWav.getDuration() + halfwaySoundWav.tStart)     && halfwaySoundWav.status === PsychoJS.Status.STARTED) {
      halfwaySoundWav.stop();  // stop the sound (if longer than duration)
      halfwaySoundWav.status = PsychoJS.Status.FINISHED;
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of halfwayBreakComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function halfwayBreakRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'halfwayBreak'-------
    for (const thisComponent of halfwayBreakComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('halfway_resp.keys', halfway_resp.keys);
    if (typeof halfway_resp.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('halfway_resp.rt', halfway_resp.rt);
        routineTimer.reset();
        }
    
    halfway_resp.stop();
    halfwayWav.stop();  // ensure sound has stopped at end of routine
    halfwaySoundWav.stop();  // ensure sound has stopped at end of routine
    // the Routine "halfwayBreak" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _valley_resp_allKeys;
var InnInstComponents;
function InnInstRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'InnInst'-------
    t = 0;
    InnInstClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    valley_resp.keys = undefined;
    valley_resp.rt = undefined;
    _valley_resp_allKeys = [];
    inn_voice.setVolume(1);
    // keep track of which components have finished
    InnInstComponents = [];
    InnInstComponents.push(valleyinst_img);
    InnInstComponents.push(valley_resp);
    InnInstComponents.push(inn_voice);
    
    for (const thisComponent of InnInstComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function InnInstRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'InnInst'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = InnInstClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *valleyinst_img* updates
    if (t >= 0.0 && valleyinst_img.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      valleyinst_img.tStart = t;  // (not accounting for frame time here)
      valleyinst_img.frameNStart = frameN;  // exact frame index
      
      valleyinst_img.setAutoDraw(true);
    }

    
    // *valley_resp* updates
    if (t >= 0.0 && valley_resp.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      valley_resp.tStart = t;  // (not accounting for frame time here)
      valley_resp.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { valley_resp.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { valley_resp.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { valley_resp.clearEvents(); });
    }

    if (valley_resp.status === PsychoJS.Status.STARTED) {
      let theseKeys = valley_resp.getKeys({keyList: ['space'], waitRelease: false});
      _valley_resp_allKeys = _valley_resp_allKeys.concat(theseKeys);
      if (_valley_resp_allKeys.length > 0) {
        valley_resp.keys = _valley_resp_allKeys[_valley_resp_allKeys.length - 1].name;  // just the last key pressed
        valley_resp.rt = _valley_resp_allKeys[_valley_resp_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // start/stop inn_voice
    if (t >= 0.0 && inn_voice.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      inn_voice.tStart = t;  // (not accounting for frame time here)
      inn_voice.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ inn_voice.play(); });  // screen flip
      inn_voice.status = PsychoJS.Status.STARTED;
    }
    if (t >= (inn_voice.getDuration() + inn_voice.tStart)     && inn_voice.status === PsychoJS.Status.STARTED) {
      inn_voice.stop();  // stop the sound (if longer than duration)
      inn_voice.status = PsychoJS.Status.FINISHED;
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of InnInstComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function InnInstRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'InnInst'-------
    for (const thisComponent of InnInstComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('valley_resp.keys', valley_resp.keys);
    if (typeof valley_resp.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('valley_resp.rt', valley_resp.rt);
        routineTimer.reset();
        }
    
    valley_resp.stop();
    inn_voice.stop();  // ensure sound has stopped at end of routine
    // the Routine "InnInst" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var charUnlockMovClock;
var charUnlockMov;
var _CharUn_resp_allKeys;
var charUnlockComponents;
function charUnlockRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'charUnlock'-------
    t = 0;
    charUnlockClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    charUnlockMovClock = new util.Clock();
    charUnlockMov = new visual.MovieStim({
      win: psychoJS.window,
      name: 'charUnlockMov',
      units: 'pix',
      movie: unlockMov,
      pos: [0, 0],
      size: [200, 200],
      ori: 0,
      opacity: 1,
      loop: true,
      noAudio: false,
      });
    CharUn_resp.keys = undefined;
    CharUn_resp.rt = undefined;
    _CharUn_resp_allKeys = [];
    CharUnlockImg.setImage(charUnlockImg);
    amountofGold.setText(coins);
    charUn_voice = new sound.Sound({
    win: psychoJS.window,
    value: charUnlockWav,
    secs: -1,
    });
    charUn_voice.setVolume(1);
    if ((blockCount === 4)) {
        continueRoutine = false;
    }
    
    // keep track of which components have finished
    charUnlockComponents = [];
    charUnlockComponents.push(charUnlockMov);
    charUnlockComponents.push(CharUn_resp);
    charUnlockComponents.push(CharUnlockImg);
    charUnlockComponents.push(amountofGold);
    charUnlockComponents.push(charUn_voice);
    
    for (const thisComponent of charUnlockComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function charUnlockRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'charUnlock'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = charUnlockClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *charUnlockMov* updates
    if (t >= 0.0 && charUnlockMov.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      charUnlockMov.tStart = t;  // (not accounting for frame time here)
      charUnlockMov.frameNStart = frameN;  // exact frame index
      
      charUnlockMov.setAutoDraw(true);
    }

    frameRemains = 0.0 + undefined - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (charUnlockMov.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      charUnlockMov.setAutoDraw(false);
    }

    
    // *CharUn_resp* updates
    if (t >= 0.0 && CharUn_resp.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      CharUn_resp.tStart = t;  // (not accounting for frame time here)
      CharUn_resp.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { CharUn_resp.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { CharUn_resp.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { CharUn_resp.clearEvents(); });
    }

    if (CharUn_resp.status === PsychoJS.Status.STARTED) {
      let theseKeys = CharUn_resp.getKeys({keyList: ['space'], waitRelease: false});
      _CharUn_resp_allKeys = _CharUn_resp_allKeys.concat(theseKeys);
      if (_CharUn_resp_allKeys.length > 0) {
        CharUn_resp.keys = _CharUn_resp_allKeys[_CharUn_resp_allKeys.length - 1].name;  // just the last key pressed
        CharUn_resp.rt = _CharUn_resp_allKeys[_CharUn_resp_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *CharUnlockImg* updates
    if (t >= 0.0 && CharUnlockImg.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      CharUnlockImg.tStart = t;  // (not accounting for frame time here)
      CharUnlockImg.frameNStart = frameN;  // exact frame index
      
      CharUnlockImg.setAutoDraw(true);
    }

    
    // *amountofGold* updates
    if (t >= 0.0 && amountofGold.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      amountofGold.tStart = t;  // (not accounting for frame time here)
      amountofGold.frameNStart = frameN;  // exact frame index
      
      amountofGold.setAutoDraw(true);
    }

    // start/stop charUn_voice
    if (t >= 0.0 && charUn_voice.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      charUn_voice.tStart = t;  // (not accounting for frame time here)
      charUn_voice.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ charUn_voice.play(); });  // screen flip
      charUn_voice.status = PsychoJS.Status.STARTED;
    }
    if (t >= (charUn_voice.getDuration() + charUn_voice.tStart)     && charUn_voice.status === PsychoJS.Status.STARTED) {
      charUn_voice.stop();  // stop the sound (if longer than duration)
      charUn_voice.status = PsychoJS.Status.FINISHED;
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of charUnlockComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function charUnlockRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'charUnlock'-------
    for (const thisComponent of charUnlockComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('CharUn_resp.keys', CharUn_resp.keys);
    if (typeof CharUn_resp.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('CharUn_resp.rt', CharUn_resp.rt);
        routineTimer.reset();
        }
    
    CharUn_resp.stop();
    charUn_voice.stop();  // ensure sound has stopped at end of routine
    // the Routine "charUnlock" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _valley_resp1_allKeys;
var valleyInstComponents;
function valleyInstRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'valleyInst'-------
    t = 0;
    valleyInstClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    valleyinst_image.setImage(valleyImg);
    valley_resp1.keys = undefined;
    valley_resp1.rt = undefined;
    _valley_resp1_allKeys = [];
    valley_voice = new sound.Sound({
    win: psychoJS.window,
    value: valleyWav,
    secs: -1,
    });
    valley_voice.setVolume(1);
    // keep track of which components have finished
    valleyInstComponents = [];
    valleyInstComponents.push(valleyinst_image);
    valleyInstComponents.push(valley_resp1);
    valleyInstComponents.push(valley_voice);
    
    for (const thisComponent of valleyInstComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function valleyInstRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'valleyInst'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = valleyInstClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *valleyinst_image* updates
    if (t >= 0.0 && valleyinst_image.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      valleyinst_image.tStart = t;  // (not accounting for frame time here)
      valleyinst_image.frameNStart = frameN;  // exact frame index
      
      valleyinst_image.setAutoDraw(true);
    }

    
    // *valley_resp1* updates
    if (t >= 0.0 && valley_resp1.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      valley_resp1.tStart = t;  // (not accounting for frame time here)
      valley_resp1.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { valley_resp1.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { valley_resp1.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { valley_resp1.clearEvents(); });
    }

    if (valley_resp1.status === PsychoJS.Status.STARTED) {
      let theseKeys = valley_resp1.getKeys({keyList: ['space'], waitRelease: false});
      _valley_resp1_allKeys = _valley_resp1_allKeys.concat(theseKeys);
      if (_valley_resp1_allKeys.length > 0) {
        valley_resp1.keys = _valley_resp1_allKeys[_valley_resp1_allKeys.length - 1].name;  // just the last key pressed
        valley_resp1.rt = _valley_resp1_allKeys[_valley_resp1_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // start/stop valley_voice
    if (t >= 0.0 && valley_voice.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      valley_voice.tStart = t;  // (not accounting for frame time here)
      valley_voice.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ valley_voice.play(); });  // screen flip
      valley_voice.status = PsychoJS.Status.STARTED;
    }
    if (t >= (valley_voice.getDuration() + valley_voice.tStart)     && valley_voice.status === PsychoJS.Status.STARTED) {
      valley_voice.stop();  // stop the sound (if longer than duration)
      valley_voice.status = PsychoJS.Status.FINISHED;
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of valleyInstComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function valleyInstRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'valleyInst'-------
    for (const thisComponent of valleyInstComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('valley_resp1.keys', valley_resp1.keys);
    if (typeof valley_resp1.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('valley_resp1.rt', valley_resp1.rt);
        routineTimer.reset();
        }
    
    valley_resp1.stop();
    valley_voice.stop();  // ensure sound has stopped at end of routine
    // the Routine "valleyInst" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _key_resp_end_allKeys;
var gameEndComponents;
function gameEndRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'gameEnd'-------
    t = 0;
    gameEndClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    end_voice.setVolume(1);
    key_resp_end.keys = undefined;
    key_resp_end.rt = undefined;
    _key_resp_end_allKeys = [];
    // keep track of which components have finished
    gameEndComponents = [];
    gameEndComponents.push(theEnd_text);
    gameEndComponents.push(endMov);
    gameEndComponents.push(EndImg);
    gameEndComponents.push(end_voice);
    gameEndComponents.push(key_resp_end);
    
    for (const thisComponent of gameEndComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function gameEndRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'gameEnd'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = gameEndClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *theEnd_text* updates
    if (t >= 0.0 && theEnd_text.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      theEnd_text.tStart = t;  // (not accounting for frame time here)
      theEnd_text.frameNStart = frameN;  // exact frame index
      
      theEnd_text.setAutoDraw(true);
    }

    
    // *endMov* updates
    if (t >= 0.0 && endMov.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      endMov.tStart = t;  // (not accounting for frame time here)
      endMov.frameNStart = frameN;  // exact frame index
      
      endMov.setAutoDraw(true);
    }

    
    // *EndImg* updates
    if (t >= 3 && EndImg.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      EndImg.tStart = t;  // (not accounting for frame time here)
      EndImg.frameNStart = frameN;  // exact frame index
      
      EndImg.setAutoDraw(true);
    }

    // start/stop end_voice
    if (t >= 0.0 && end_voice.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      end_voice.tStart = t;  // (not accounting for frame time here)
      end_voice.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ end_voice.play(); });  // screen flip
      end_voice.status = PsychoJS.Status.STARTED;
    }
    if (t >= (end_voice.getDuration() + end_voice.tStart)     && end_voice.status === PsychoJS.Status.STARTED) {
      end_voice.stop();  // stop the sound (if longer than duration)
      end_voice.status = PsychoJS.Status.FINISHED;
    }
    
    // *key_resp_end* updates
    if (t >= 0.0 && key_resp_end.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      key_resp_end.tStart = t;  // (not accounting for frame time here)
      key_resp_end.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { key_resp_end.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { key_resp_end.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { key_resp_end.clearEvents(); });
    }

    if (key_resp_end.status === PsychoJS.Status.STARTED) {
      let theseKeys = key_resp_end.getKeys({keyList: ['space'], waitRelease: false});
      _key_resp_end_allKeys = _key_resp_end_allKeys.concat(theseKeys);
      if (_key_resp_end_allKeys.length > 0) {
        key_resp_end.keys = _key_resp_end_allKeys[_key_resp_end_allKeys.length - 1].name;  // just the last key pressed
        key_resp_end.rt = _key_resp_end_allKeys[_key_resp_end_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of gameEndComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function gameEndRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'gameEnd'-------
    for (const thisComponent of gameEndComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    end_voice.stop();  // ensure sound has stopped at end of routine
    psychoJS.experiment.addData('key_resp_end.keys', key_resp_end.keys);
    if (typeof key_resp_end.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('key_resp_end.rt', key_resp_end.rt);
        routineTimer.reset();
        }
    
    key_resp_end.stop();
    // the Routine "gameEnd" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _credits_resp_allKeys;
var CreditsComponents;
function CreditsRoutineBegin(trials) {
  return function () {
    //------Prepare to start Routine 'Credits'-------
    t = 0;
    CreditsClock.reset(); // clock
    frameN = -1;
    // update component parameters for each repeat
    credits_wav.setVolume(1);
    credits_resp.keys = undefined;
    credits_resp.rt = undefined;
    _credits_resp_allKeys = [];
    // keep track of which components have finished
    CreditsComponents = [];
    CreditsComponents.push(Credits_text);
    CreditsComponents.push(credits_wav);
    CreditsComponents.push(credits_resp);
    
    for (const thisComponent of CreditsComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    
    return Scheduler.Event.NEXT;
  };
}


function CreditsRoutineEachFrame(trials) {
  return function () {
    //------Loop for each frame of Routine 'Credits'-------
    let continueRoutine = true; // until we're told otherwise
    // get current time
    t = CreditsClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *Credits_text* updates
    if (t >= 0.0 && Credits_text.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      Credits_text.tStart = t;  // (not accounting for frame time here)
      Credits_text.frameNStart = frameN;  // exact frame index
      
      Credits_text.setAutoDraw(true);
    }

    // start/stop credits_wav
    if (t >= 0.0 && credits_wav.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      credits_wav.tStart = t;  // (not accounting for frame time here)
      credits_wav.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ credits_wav.play(); });  // screen flip
      credits_wav.status = PsychoJS.Status.STARTED;
    }
    if (t >= (credits_wav.getDuration() + credits_wav.tStart)     && credits_wav.status === PsychoJS.Status.STARTED) {
      credits_wav.stop();  // stop the sound (if longer than duration)
      credits_wav.status = PsychoJS.Status.FINISHED;
    }
    
    // *credits_resp* updates
    if (t >= 0.0 && credits_resp.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      credits_resp.tStart = t;  // (not accounting for frame time here)
      credits_resp.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { credits_resp.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { credits_resp.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { credits_resp.clearEvents(); });
    }

    if (credits_resp.status === PsychoJS.Status.STARTED) {
      let theseKeys = credits_resp.getKeys({keyList: ['space'], waitRelease: false});
      _credits_resp_allKeys = _credits_resp_allKeys.concat(theseKeys);
      if (_credits_resp_allKeys.length > 0) {
        credits_resp.keys = _credits_resp_allKeys[_credits_resp_allKeys.length - 1].name;  // just the last key pressed
        credits_resp.rt = _credits_resp_allKeys[_credits_resp_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of CreditsComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function CreditsRoutineEnd(trials) {
  return function () {
    //------Ending Routine 'Credits'-------
    for (const thisComponent of CreditsComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    credits_wav.stop();  // ensure sound has stopped at end of routine
    psychoJS.experiment.addData('credits_resp.keys', credits_resp.keys);
    if (typeof credits_resp.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('credits_resp.rt', credits_resp.rt);
        routineTimer.reset();
        }
    
    credits_resp.stop();
    // the Routine "Credits" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


function endLoopIteration(thisScheduler, loop) {
  // ------Prepare for next entry------
  return function () {
    if (typeof loop !== 'undefined') {
      // ------Check if user ended loop early------
      if (loop.finished) {
        // Check for and save orphaned data
        if (psychoJS.experiment.isEntryEmpty()) {
          psychoJS.experiment.nextEntry(loop);
        }
      thisScheduler.stop();
      } else {
        const thisTrial = loop.getCurrentTrial();
        if (typeof thisTrial === 'undefined' || !('isTrials' in thisTrial) || thisTrial.isTrials) {
          psychoJS.experiment.nextEntry(loop);
        }
      }
    return Scheduler.Event.NEXT;
    }
  };
}


function importConditions(trials) {
  return function () {
    psychoJS.importAttributes(trials.getCurrentTrial());
    return Scheduler.Event.NEXT;
    };
}


function quitPsychoJS(message, isCompleted) {
  // Check for and save orphaned data
  if (psychoJS.experiment.isEntryEmpty()) {
    psychoJS.experiment.nextEntry();
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  psychoJS.window.close();
  psychoJS.quit({message: message, isCompleted: isCompleted});
  
  return Scheduler.Event.QUIT;
}
