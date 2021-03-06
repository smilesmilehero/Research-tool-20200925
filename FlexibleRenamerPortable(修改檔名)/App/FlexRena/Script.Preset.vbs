Public Const vbNone		= &h0000 '//No flags.

'//Arrange dight of number
Public Const vbSignWithinWidth	= &h0001 '//Includes sign whthin the given Width.
Public Const vbSignAlways		= &h0002 '//Prefix the output value with a sign(+ or -).
Function DigitNumber(Number, Width, SignOption)
	DigitNumber = CStr(CLng(Abs(Number)))
	IsSignWithinDigit = CBool(SignOption And vbSignWithinWidth)
	IsSignAlways = CBool(SignOption And vbSignAlways)
	IsNegativeSigned = CBool(Number<0)
	If (IsSignWithinDigit And (IsSignAlways Or IsNegativeSigned)) Then Width = Width-1
	nZero = Width-Len(DigitNumber)
	If (nZero>0) Then DigitNumber = String(nZero, "0") & DigitNumber
	If (IsNegativeSigned) Then
		DigitNumber = "-" & DigitNumber
	ElseIf (IsSignAlways) Then
		DigitNumber = "+" & DigitNumber
	End If
End Function

'//Halfwidth <-> Fullwidth
FullKana = Split("ヴ ガ ギ グ ゲ ゴ ザ ジ ズ ゼ ゾ ダ ヂ ヅ デ ド バ ビ ブ ベ ボ パ ピ プ ペ ポ ゛ 。 「 」 、 ・ ヲ ァ ィ ゥ ェ ォ ャ ュ ョ ッ ー ア イ ウ エ オ カ キ ク ケ コ サ シ ス セ ソ タ チ ツ テ ト ナ ニ ヌ ネ ノ ハ ヒ フ ヘ ホ マ ミ ム メ モ ヤ ユ ヨ ラ リ ル レ ロ ワ ン ゛ ゜")
HalfKana = Split("ｳﾞ ｶﾞ ｷﾞ ｸﾞ ｹﾞ ｺﾞ ｻﾞ ｼﾞ ｽﾞ ｾﾞ ｿﾞ ﾀﾞ ﾁﾞ ﾂﾞ ﾃﾞ ﾄﾞ ﾊﾞ ﾋﾞ ﾌﾞ ﾍﾞ ﾎﾞ ﾊﾟ ﾋﾟ ﾌﾟ ﾍﾟ ﾎﾟ ﾞ ｡ ｢ ｣ ､ ･ ｦ ｧ ｨ ｩ ｪ ｫ ｬ ｭ ｮ ｯ ｰ ｱ ｲ ｳ ｴ ｵ ｶ ｷ ｸ ｹ ｺ ｻ ｼ ｽ ｾ ｿ ﾀ ﾁ ﾂ ﾃ ﾄ ﾅ ﾆ ﾇ ﾈ ﾉ ﾊ ﾋ ﾌ ﾍ ﾎ ﾏ ﾐ ﾑ ﾒ ﾓ ﾔ ﾕ ﾖ ﾗ ﾘ ﾙ ﾚ ﾛ ﾜ ﾝ ﾞ ﾟ")
FullAscii = "！”＃＄％＆’（）＊＋，－．／０１２３４５６７８９：；＜＝＞？＠ＡＢＣＤＥＦＧＨＩＪＫＬＭＮＯＰＱＲＳＴＵＶＷＸＹＺ［￥］ａｂｃｄｅｆｇｈｉｊｋｌｍｎｏｐｑｒｓｔｕｖｗｘｙｚ｛｜｝～　"
HalfAscii = "!""#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]abcdefghijklmnopqrstuvwxyz{|}~ "
Function toFullKana(str)
	toFullKana = str
	For i = LBound(HalfKana) to UBound(HalfKana)	: toFullKana = Replace(toFullKana, HalfKana(i), FullKana(i),1,-1,vbBinaryCompare)	: Next
End Function
Function toHalfKana(str)
	toHalfKana = str
	For i = LBound(FullKana) to UBound(FullKana)	: toHalfKana = Replace(toHalfKana, FullKana(i), HalfKana(i),1,-1,vbBinaryCompare)	: Next
End Function
Function toFull(str)
	toFull = toFullKana(str)
	For i = 1 to Len(HalfAscii)	: toFull = Replace(toFull, Mid(HalfAscii, i, 1), Mid(FullAscii, i, 1),1,-1,vbBinaryCompare)	: Next
End Function
Function toHalf(str)
	toHalf = toHalfKana(str)
	For i = 1 to Len(FullAscii)	: toHalf = Replace(toHalf , Mid(FullAscii, i, 1), Mid(HalfAscii, i, 1),1,-1,vbBinaryCompare)	: Next
End Function

'//RomanKana
HiraGana = Split("ヴ が ぎ ぐ げ ご ざ じ ず ぜ ぞ だ ぢ づ で ど ば び ぶ べ ぼ ぱ ぴ ぷ ぺ ぽ ゛ 。 「 」 、 ・ を ぁ ぃ ぅ ぇ ぉ ゃ ゅ ょ っ ー あ い う え お か き く け こ さ し す せ そ た ち つ て と な に ぬ ね の は ひ ふ へ ほ ま み む め も や ゆ よ ら り る れ ろ わ ん ゛ ゜")
'------------------------------------------------------------------------------------------------------------------------------------------
'このローマ字表は、ヘボン式・英国式・その他の表記を参考に、渡辺真 氏が編纂したものです。(http://makotowatana.ld.infoseek.co.jp/hotvbs.html)
KataKana = Split("ッビャー ッビェー ッビョー ッビュー ッディー ッドゥー ッデュー ッグァー ッグェー ッグィー ッグォー ッギャー ッギェー ッギョー ッギュー ッジャー ッジェー ッジョー ッジュー ッンバー ッンベー ッンビー ッンボー ッンブー ッンパー ッンペー ッンピー ッンポー ッンプー ッピャー ッピェー ッピョー ッピュー ッヴァー ッヴェー ッヴィー ッヴォー ッヴョー ッヴュー ッズィー ッヂャー ッヂョー ッヂュー ッファー ッフェー ッフィー ッフォー ッフャー ッフョー ッフュー ッヒャー ッヒェー ッヒョー ッヒュー ックァー ックェー ックィー ックォー ッキャー ッキェー ッキョー ッキュー ッンマー ッンメー ッンミー ッンモー ッンムー ッミャー ッミェー ッミョー ッミュー ッンアー ッンエー ッンイー ッンナー ッンネー ッンニー ッンノー ッンヌー ッンオー ッンウー ッニャー ッニェー ッニェー ッニョー ッニュー ッリャー ッリェー ッリョー ッリュー ッシャー ッシェー ッシェー ッショー ッシュー ッスィー ッチャー ッチェー ッチェー ッチョー ッチュー ッティー ッオオー ッオウー ッツァー ッツェー ッツィー ッツォー ットゥー ッテュー ッウェー ッウィー ッウォー ッイェー ッバー ッベー ッビー ッボー ッブー ッビャ ッビェ ッビョ ッビュ ビャー ビェー ビョー ビュー ッダー ッデー ッディ ッドー ッドゥ ッデュ ッヅー ディー ドゥー デュー ッガー ッゲー ッギー ッゴー ッグー ッグァ ッグェ ッグィ ッグォ ッギャ ッギェ ッギョ ッギュ グァー グェー グィー グォー ギャー ギェー ギョー ギュー ジャー ジェー ッジャ ッジェ ッジー ッジョ ッジュ ジョー ジュー ンバー ンベー ンビー ンボー ンブー ッンバ ッンベ ッンビ ッンボ ッンブ ッンパ ッンペ ッンピ ッンポ ッンプ ンパー ンペー ンピー ンポー ンプー ッパー ッペー ッピー ッポー ップー ッピャ ッピェ ッピョ ッピュ ピャー ピェー ピョー ピュー ヴァー ヴェー ヴィー ヴォー ッヴァ ッヴェ ッヴィ ッヴォ ッヴー ッヴョ ッヴュ ヴョー ヴュー ズィー ヂャー ヂョー ヂュー ッザー ッゼー ッズィ ッヂー ッゾー ッズー ッヂャ ッヂョ ッヂュ チャー チェー チェー チョー チュー ファー フェー ッファ ッフェ ッフィ ッフォ ッフー ッフャ ッフョ ッフュ フィー フォー フャー フョー フュー ッハー ッヘー ッヒー ッホー ッヒャ ッヒェ ッヒョ ッヒュ ヒャー ヒェー ヒョー ヒュー ッカー ッケー ッキー ッコー ックー ックァ ックェ ックィ ックォ ッキャ ッキェ ッキョ ッキュ クァー クェー クィー クォー キャー キェー キョー キュー ンマー ッマー ンメー ッメー ンミー ッミー ッンマ ッンメ ッンミ ッンモ ッンム ンモー ッモー ンムー ッムー ッミャ ッミェ ッミョ ッミュ ミャー ミェー ミョー ミュー ンアー ッナー ンエー ッネー ンイー ッニー ッンア ンナー ッンエ ンネー ッンイ ンニー ッンナ ッンネ ッンニ ッンノ ッンヌ ッンオ ンノー ッンウ ンヌー ンオー ッノー ンウー ッヌー ッニャ ッニェ ッニョ ッニュ ニャー ニェー ニョー ニュー オオー オウー ッラー ッレー ッリー ッロー ッルー ッリャ ッリェ ッリョ ッリュ リャー リェー リョー リュー シャー シェー ショー シュー スィー ッサー ッセー ッシャ ッシェ ッシー ッショ ッシュ ッスィ ッソー ッスー ッチャ ッチェ ッチェ ッチー ッチョ ッチュ ティー ツァー ツェー ツィー ツォー ッアー ッター ッエー ッテー ッティ ッイー ッオオ ッオウ ッオー ットー ッツァ ッツェ ッツィ ッツォ ッツー ットゥ ッウー ッテュ トゥー テュー ウェー ウィー ウォー ッワー ッウェ ッヱー ッウィ ッヰー ッウォ ッヲー イェー イェー ッヤー ッイェ ッイェ ッヨー ッユー バー ッバ ッベ ッビ ッボ ッブ ベー ビー ボー ブー ビャ ビェ ビィ ビョ ビュ ダー ッダ ッデ ッド ッヅ デー ディ ドー ドゥ ドォ デャ デョ デュ ヅー ガー ゲー ッガ ッゲ ッギ ッゴ ッグ ギー ゴー グー グァ グェ グィ グォ ギャ ギェ ギィ ギョ ギュ ジャ ジェ ジー ッジ ジョ ジュ ジィ ンバ ンベ ンビ ンボ ンブ ンパ ンペ ンピ ンポ ンプ パー ペー ピー ポー ッパ ッペ ッピ ッポ ップ プー ピャ ピェ ピィ ピョ ピュ ヴァ ヴェ ヴィ ヴォ ヴー ッヴ ヴャ ヴョ ヴュ ザー ゼー ズィ ヂー ゾー ズー ヂャ ヂョ ヂュ ッザ ッゼ ッヂ ッゾ ッズ アー チャ チェ チー チョ チュ エー ファ フェ ッフ フィ フォ フー フゥ フャ フョ フュ ハー ヘー ッハ ッヘ ッヒ ッホ ヒー ホー ヒャ ヒェ ヒィ ヒョ ヒュ イー カー ケー キー ッカ ッケ ッキ ッコ ック コー クー クァ クェ クィ クォ キャ キェ キィ キョ キュ マー メー ミー ンマ ッマ ンメ ッメ ンミ ッミ ンモ ッモ ンム ッム モー ムー ミャ ミェ ミィ ミョ ミュ ナー ネー ニー ンア ッナ ンエ ッネ ンイ ッニ ンナ ンネ ンニ ンノ ンヌ ンオ ッノ ンウ ッヌ ノー ヌー ニャ ニェ ニィ ニョ ニュ オオ オウ オー クャ クョ クュ ラー レー リー ロー ッラ ッレ ッリ ッロ ッル ルー リャ リェ リィ リョ リュ サー セー シャ シェ シー ショ シュ スィ ソー ッサ ッセ ッシ ッソ ッス スー スァ スォ スゥ ター ッチ テー ティ トー ツァ ツェ ツィ ツォ ツー ッア ッタ ッエ ッテ ッイ ッオ ット ッツ ッウ トゥ テョ テュ ウー ワー ウェ ヱー ウィ ヰー ウォ ヲー ッワ ッヱ ッヰ ッヲ ヤー イェ ヨー ユー ッヤ ッヨ ッユ バ ベ ビ ボ ブ ダ デ ド ヅ ガ ゲ ギ ゴ グ ジ パ ペ ピ ポ プ ヴ ザ ゼ ヂ ゾ ズ ー ア ァ チ エ ェ フ ハ ヘ ヒ ホ イ ィ カ ケ キ コ ク マ メ ミ モ ム ン ナ ネ ニ ノ ヌ オ ォ ラ レ リ ロ ル サ セ シ ソ ス ッ タ テ ト ツ ウ ゥ ワ ヱ ヰ ヲ ヤ ャ ヨ ョ ユ ュ")
RomanKana = Split("bbyaa bbyee bbyoo bbyuu ddii dduu ddyuu ggwaa ggwee ggwii ggwoo ggyaa ggyee ggyoo ggyuu jjaa jjee jjoo jjuu mmbaa mmbee mmbii mmboo mmbuu mmpaa mmpee mmpii mmpoo mmpuu ppyaa ppyee ppyoo ppyuu vvaa vvee vvii vvoo vvyoo vvyuu zz'ii zzyaa zzyoo zzyuu ffaa ffee ffii ffoo ffyaa ffyoo ffyuu hhyaa hhyee hhyoo hhyuu kkwaa kkwee kkwii kkwoo kkyaa kkyee kkyoo kkyuu mmmaa mmmee mmmii mmmoo mmmuu mmyaa mmyee mmyoo mmyuu nnnaa nnnee nnnii nnnnaa nnnnee nnnnii nnnnoo nnnnuu nnnoo nnnuu nnyaa nnyee nnyee nnyoo nnyuu rryaa rryee rryoo rryuu sshaa sshee sshee sshoo sshuu ss'ii tchaa tchee tchee tchoo tchuu ttii ttohh ttohh ttsaa ttsee ttsii ttsoo ttuu ttyuu wwee wwii wwoo yyee bbaa bbee bbii bboo bbuu bbya bbye bbyo bbyu byaa byee byoo byuu ddaa ddee ddi ddoo ddu ddyu ddzuu dii duu dyuu ggaa ggee ggii ggoo gguu ggwa ggwe ggwi ggwo ggya ggye ggyo ggyu gwaa gwee gwii gwoo gyaa gyee gyoo gyuu jaa jee jja jje jjii jjo jju joo juu mbaa mbee mbii mboo mbuu mmba mmbe mmbi mmbo mmbu mmpa mmpe mmpi mmpo mmpu mpaa mpee mpii mpoo mpuu ppaa ppee ppii ppoo ppuu ppya ppye ppyo ppyu pyaa pyee pyoo pyuu vaa vee vii voo vva vve vvi vvo vvuu vvyo vvyu vyoo vyuu z'ii zyaa zyoo zyuu zzaa zzee zz'i zzjii zzoo zzuu zzya zzyo zzyu chaa chee chee choo chuu faa fee ffa ffe ffi ffo ffuu ffya ffyo ffyu fii foo fyaa fyoo fyuu hhaa hhee hhii hhoo hhya hhye hhyo hhyu hyaa hyee hyoo hyuu kkaa kkee kkii kkoo kkuu kkwa kkwe kkwi kkwo kkya kkye kkyo kkyu kwaa kwee kwii kwoo kyaa kyee kyoo kyuu mmaa mmaa mmee mmee mmii mmii mmma mmme mmmi mmmo mmmu mmoo mmoo mmuu mmuu mmya mmye mmyo mmyu myaa myee myoo myuu nnaa nnaa nnee nnee nnii nnii nnna nnnaa nnne nnnee nnni nnnii nnnna nnnne nnnni nnnno nnnnu nnno nnnoo nnnu nnnuu nnoo nnoo nnuu nnuu nnya nnye nnyo nnyu nyaa nyee nyoo nyuu ohh ohh rraa rree rrii rroo rruu rrya rrye rryo rryu ryaa ryee ryoo ryuu shaa shee shoo shuu s'ii ssaa ssee ssha sshe sshii ssho sshu ss'i ssoo ssuu tcha tche tche tchii tcho tchu tii tsaa tsee tsii tsoo ttaa ttaa ttee ttee tti ttii ttoh ttoh ttoo ttoo ttsa ttse ttsi ttso ttsuu ttu ttuu ttyu tuu tyuu wee wii woo wwaa wwe wwee wwi wwii wwo wwoo yee yee yyaa yye yye yyoo yyuu baa bba bbe bbi bbo bbu bee bii boo buu bya bye byi byo byu daa dda dde ddo ddzu dee di doo du dwo dya dyo dyu dzuu gaa gee gga gge ggi ggo ggu gii goo guu gwa gwe gwi gwo gya gye gyi gyo gyu ja je jii jji jo ju jyi mba mbe mbi mbo mbu mpa mpe mpi mpo mpu paa pee pii poo ppa ppe ppi ppo ppu puu pya pye pyi pyo pyu va ve vi vo vuu vvu vya vyo vyu zaa zee z'i zjii zoo zuu zya zyo zyu zza zze zzji zzo zzu aa cha che chii cho chu ee fa fe ffu fi fo fuu fwu fya fyo fyu haa hee hha hhe hhi hho hii hoo hya hye hyi hyo hyu ii kaa kee kii kka kke kki kko kku koo kuu kwa kwe kwi kwo kya kye kyi kyo kyu maa mee mii mma mma mme mme mmi mmi mmo mmo mmu mmu moo muu mya mye myi myo myu naa nee nii nna nna nne nne nni nni nnna nnne nnni nnno nnnu nno nno nnu nnu noo nuu nya nye nyi nyo nyu oh oh oo qya qyo qyu raa ree rii roo rra rre rri rro rru ruu rya rye ryi ryo ryu saa see sha she shii sho shu s'i soo ssa sse sshi sso ssu suu swa swo swu taa tchi tee ti too tsa tse tsi tso tsuu tta tta tte tte tti tto tto ttsu ttu tu tyo tyu uu waa we wee wi wii wo woo wwa wwe wwi wwo yaa ye yoo yuu yya yyo yyu ba be bi bo bu da de do dzu ga ge gi go gu ji pa pe pi po pu vu za ze zi zo zu - a a chi e e fu ha he hi ho i i ka ke ki ko ku ma me mi mo mu n na ne ni no nu o o ra re ri ro ru sa se shi so su t ta te to tsu u u wa we wi wo ya ya yo yo yu yu")
'------------------------------------------------------------------------------------------------------------------------------------------
Function HiraganatoKatakana(str)
	HiraganatoKatakana = str
	For i = LBound(HiraGana) to UBound(HiraGana)	: HiraganatoKatakana = Replace(HiraganatoKatakana, HiraGana(i), FullKana(i),1,-1,vbBinaryCompare)	: Next
End Function
Function toRomanKana(str)
	toRomanKana = toFullKana(HiraganatoKatakana(str))
	For i = LBound(KataKana) to UBound(KataKana)	: toRomanKana = Replace(toRomanKana, KataKana(i), RomanKana(i),1,-1,vbBinaryCompare)	: Next
End Function

'// for StrFunc.dll
Public Const vbUpperCase		= &h0001 '//Converts the string to uppercase characters.
Public Const vbLowerCase		= &h0002 '//Converts the string to lowercase characters.
Public Const vbProperCase		= &h0003 '//Converts the first letter of every word in string to uppercase.
Public Const vbWide		= &h0004 '//Converts narrow (half-width) characters in the string to wide (full-width) characters.
Public Const vbNarrow		= &h0008 '//Converts wide (full-width) characters in the string to narrow (half-width) characters.
Public Const vbKatakana		= &h0010 '//Converts Hiragana characters in the string to Katakana characters.
Public Const vbHiragana		= &h0020 '//Converts Katakana characters in the string to Hiragana characters.
Public Const vbSimplifiedChinese	= &h0100 '//Converts Traditional Chinese characters to Simplified Chinese.
Public Const vbTraditionalChinese	= &h0200 '//Converts Simplified Chinese characters to Traditional Chinese.
Public Const vbLinguisticCasing	= &h0400 '//Uses linguistic rules for casing, rather than File System (default). Valid with vbUpperCase and vbLowerCase only.
Public Const LCMAP_LOWERCASE	= &h00000100  '//Map all characters to lowercase.
Public Const LCMAP_UPPERCASE	= &h00000200  '//Map all characters to uppercase.
Public Const LCMAP_PROPERCASE	= &h00000300  '//Map the first letter of every word in string to uppercase.
Public Const LCMAP_SORTKEY	= &h00000400  '//Produce a normalized wide character sort key.
Public Const LCMAP_BYTEREV	= &h00000800  '//Use byte reversal.
Public Const LCMAP_HIRAGANA	= &h00100000  '//Map all Katakana characters to Hiragana.
Public Const LCMAP_KATAKANA	= &h00200000  '//Map all Hiragana characters to Katakana.
Public Const LCMAP_HALFWIDTH	= &h00400000  '//Use narrow characters where applicable.
Public Const LCMAP_FULLWIDTH	= &h00800000  '//Use wide characters where applicable.
Public Const LCMAP_LINGUISTIC_CASING	= &h01000000  '//Use linguistic rules for casing, instead of file system rules (default).
Public Const LCMAP_SIMPLIFIED_CHINESE	= &h02000000  '//Map traditional Chinese characters to simplified Chinese characters.
Public Const LCMAP_TRADITIONAL_CHINESE	= &h04000000  '//Map simplified Chinese characters to traditional Chinese characters.
Function StrFunc()			: Set StrFunc = CreateObject("StrFunc.Wrapper")						: End Function
Function StrConv(SrcStr, Conversion)	: StrConv = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, Conversion)			: End Function
Function toLowerCase(SrcStr)		: toLowerCase = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbLowerCase)		: End Function
Function toUpperCase(SrcStr)		: toUpperCase = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbUpperCase)		: End Function
Function toProperCase(SrcStr)	: toProperCase = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbProperCase)		: End Function
Function toWide(SrcStr)		: toWide = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbWide)			: End Function
Function toNarrow(SrcStr)		: toNarrow = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbNarrow)			: End Function
Function toKatakana(SrcStr)		: toKatakana = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbKatakana)			: End Function
Function toHiragana(SrcStr)		: toHiragana = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbHiragana)			: End Function
Function toSimplifiedChinese(SrcStr)	: toSimplifiedChinese = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbSimplifiedChinese)	: End Function
Function toTraditionalChinese(SrcStr)	: toTraditionalChinese = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbTraditionalChinese)	: End Function
Function toLinguisticCasing(SrcStr)	: toLinguisticCasing = CreateObject("StrFunc.Wrapper").StrConv(SrcStr, vbLinguisticCasing)	: End Function
