基于OpenFOAM-9开发，主求解器为neutronicThermionicFoam.C，
运行时需要搭配OpenFOAM自带常微分方程求解方法：
默认中子计算参数，包括缓发中子份额、中子代时间、衰变常数等，存储在neutronODE.H中；
如需引入电特性计算，则需要为发射极下的geometricFParameter中electrod关键词命名为>1的数，接收极下geometricFParameter中electrod关键词命名为<-1的数，每根TFE对应一组电极。
目前可选择固定输出电流/固定输出电压/固定输出电阻计算模式，在createTFEFields.H中进行修改，


运行示例：
主求解器为neutronicThermionicFoam rodas34

目前仅支持CPU单核运行
